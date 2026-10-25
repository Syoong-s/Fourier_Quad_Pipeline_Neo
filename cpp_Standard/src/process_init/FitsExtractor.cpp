#include "process_init/FitsExtractor.hpp"
#include "InitConfig.hpp"

#include <fitsio.h>
#include <algorithm>
#include <array>
#include <cerrno>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fqinit {
namespace {
namespace fs = std::filesystem;
constexpr std::size_t kCopyBufferBytes = 16U * 1024U * 1024U;
constexpr std::uint64_t kFitsBlock = 2880;
constexpr unsigned kMaxHeaderBlocks = 8192;

struct PlannedImage {
    int hdu = 0;
    int ccdnum = 0;
    fs::path final_path;
    fs::path staged_path;
    bool isolated = false;
    bool sanitize = false;
    std::uint64_t offset = 0;
    std::uint64_t header_length = 0;
    std::uint64_t total_length = 0;
};
struct RawHdu {
    int hdu = 0;
    std::uint64_t offset = 0, header_length = 0, total_length = 0;
    bool image = false, compressed = false, unsafe = false, sanitize = false;
    int naxis = 0;
    int ccdnum = 0;
    bool has_ccdnum = false;
};

std::string fitsDiagnostic(int status) {
    char msg[FLEN_STATUS] = {};
    fits_get_errstatus(status, msg);
    std::ostringstream out;
    out << msg;
    char detail[FLEN_ERRMSG] = {};
    while (fits_read_errmsg(detail)) out << "; " << detail;
    return out.str();
}
void closeFits(fitsfile*& f) {
    if (!f) return;
    int status = 0;
    fits_close_file(f, &status);
    f = nullptr;
}
void replaceAll(std::string& s, const std::string& from, const std::string& to) {
    if (from.empty()) return;
    std::size_t p = 0;
    while ((p = s.find(from, p)) != std::string::npos) {
        s.replace(p, from.size(), to);
        p += to.size();
    }
}
std::string trim(const std::string& s) {
    auto b = s.find_first_not_of(' ');
    if (b == std::string::npos) return {};
    return s.substr(b, s.find_last_not_of(' ') - b + 1);
}
std::pair<int, std::size_t> imageDatatype(int bitpix) {
    switch (bitpix) {
        case BYTE_IMG: return {TBYTE, 1};
        case SBYTE_IMG: return {TSBYTE, 1};
        case SHORT_IMG: return {TSHORT, 2};
        case USHORT_IMG: return {TUSHORT, 2};
        case LONG_IMG: return {TINT, 4};
        case ULONG_IMG: return {TUINT, 4};
        case LONGLONG_IMG: return {TLONGLONG, 8};
        case ULONGLONG_IMG: return {TULONGLONG, 8};
        case FLOAT_IMG: return {TFLOAT, 4};
        case DOUBLE_IMG: return {TDOUBLE, 8};
        default: throw std::runtime_error("unsupported BITPIX " + std::to_string(bitpix));
    }
}

// Reject only provably malformed OPTIONAL card syntax. For custom keywords,
// preserve any valid FITS representation; do not try to interpret their meaning.
bool isKeyword(const std::string& key) {
    if (key.empty() || key.size() > 8) return false;
    for (unsigned char c : key) {
        if (!((c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') ||
              c == '_' || c == '-')) return false;
    }
    return true;
}
bool validNumber(std::string s) {
    s = trim(s);
    if (s.empty()) return false;
    std::size_t i = 0, nd = 0;
    if (s[i] == '+' || s[i] == '-') ++i;
    while (i < s.size() && s[i] >= '0' && s[i] <= '9') { ++i; ++nd; }
    if (i < s.size() && s[i] == '.') {
        ++i;
        while (i < s.size() && s[i] >= '0' && s[i] <= '9') { ++i; ++nd; }
    }
    if (!nd) return false;
    if (i < s.size() && (s[i] == 'E' || s[i] == 'D' || s[i] == 'e' || s[i] == 'd')) {
        ++i;
        if (i < s.size() && (s[i] == '+' || s[i] == '-')) ++i;
        std::size_t start = i;
        while (i < s.size() && s[i] >= '0' && s[i] <= '9') ++i;
        if (start == i) return false;
    }
    return i == s.size();
}
bool validInteger(std::string s) {
    s = trim(s);
    if (s.empty()) return false;
    std::size_t i = 0;
    if (s[i] == '+' || s[i] == '-') ++i;
    if (i == s.size()) return false;
    for (; i < s.size(); ++i) if (s[i] < '0' || s[i] > '9') return false;
    return true;
}
bool critical(const std::string& k) {
    static const std::set<std::string> names = {
        "SIMPLE", "XTENSION", "BITPIX", "NAXIS", "PCOUNT", "GCOUNT",
        "TFIELDS", "THEAP", "ZIMAGE", "ZCMPTYPE", "ZBITPIX", "ZNAXIS",
        "ZQUANTIZ", "ZDITHER0", "CCDNUM", "BLANK", "BSCALE", "BZERO",
        "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CD1_1", "CD1_2", "CD2_1", "CD2_2"
    };
    if (names.count(k)) return true;
    for (const char* p : {"NAXIS", "ZNAXIS", "ZTILE", "ZNAME", "ZVAL", "TTYPE",
                          "TFORM", "TSCAL", "TZERO", "TNULL"})
        if (k.rfind(p, 0) == 0) return true;
    return false;
}
bool invalidOptionalCard(const std::string& card) {
    if (card.size() != 80) return true;
    for (unsigned char c : card) if (c < 32 || c > 126) return true;
    const std::string key = trim(card.substr(0, 8));
    if (key.empty() || key == "COMMENT" || key == "HISTORY" ||
        key == "CONTINUE" || key == "HIERARCH" || key == "END") return false;
    if (!isKeyword(key) || card[8] != '=' || card[9] != ' ') return true;
    const std::string value = trim(card.substr(10));
    if (value.empty()) return false;  // FITS permits an undefined value.
    if (value.front() == '\'') {
        for (std::size_t i = card.find('\'', 10) + 1; i < card.size(); ++i) {
            if (card[i] != '\'') continue;
            if (i + 1 < card.size() && card[i+1] == '\'') { ++i; continue; }
            return false;
        }
        return true;
    }
    const std::string token = trim(value.substr(0, value.find('/')));
    if (token == "T" || token == "F") return false;
    if (token.size() >= 5 && token.front() == '(' && token.back() == ')') {
        auto comma = token.find(',');
        return comma == std::string::npos ||
               !validNumber(token.substr(1, comma-1)) ||
               !validNumber(token.substr(comma+1, token.size()-comma-2));
    }
    return !validNumber(token);
}

// Native FITS type, not CFITSIO's permissive type conversions: '7' is not an
// integer-valued FITS card, despite fits_read_key(...TINT) accepting it.
bool nativeKey(fitsfile* f, const char* key, bool integer, bool optional = false) {
    int s = 0;
    char card[FLEN_CARD] = {};
    fits_read_card(f, key, card, &s);
    if (s == KEY_NO_EXIST && optional) { fits_clear_errmsg(); return true; }
    if (s != 0) { fitsDiagnostic(s); return false; }
    const std::string c(card, std::strlen(card));
    if (c.size() < 11 || c[8] != '=' || c[9] != ' ') return false;
    const auto raw = trim(c.substr(10));
    if (raw.empty() || raw.front() == '\'') return false;
    const auto token = trim(raw.substr(0, raw.find('/')));
    return integer ? validInteger(token) : validNumber(token);
}
void validateRequired(fitsfile* f, ProductKind kind) {
    int status = 0, naxis = 0, ccdnum = 0;
    fits_get_img_dim(f, &naxis, &status);
    if (status || naxis != 2) throw std::runtime_error("required NAXIS != 2");
    if (!nativeKey(f, InitConfig::CCDNUM_KEYWORD, true))
        throw std::runtime_error("required CCDNUM is missing or not a FITS integer");
    status = 0;
    fits_read_key(f, TINT, InitConfig::CCDNUM_KEYWORD, &ccdnum, nullptr, &status);
    if (status) throw std::runtime_error("CCDNUM cannot be represented as int: " + fitsDiagnostic(status));
    if (kind == ProductKind::Science) {
        for (const char* k : {"CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2",
                              "CD1_1", "CD1_2", "CD2_1", "CD2_2"})
            if (!nativeKey(f, k, false))
                throw std::runtime_error(std::string("required WCS has wrong FITS type: ") + k);
    }
    for (const char* k : {"BSCALE", "BZERO"})
        if (!nativeKey(f, k, false, true))
            throw std::runtime_error(std::string("invalid pixel scaling type: ") + k);
    if (!nativeKey(f, "BLANK", true, true))
        throw std::runtime_error("invalid BLANK type");
}

struct HeaderFailure : std::runtime_error { using std::runtime_error::runtime_error; };
// Same normal extraction as Neo: no pixel prescan, 16MiB read/write chunks.
void copyCurrentImage(fitsfile* input, const fs::path& output_path, bool relaxed,
                      const fs::path& source, int physical_hdu) {
    int status = 0, bitpix = 0, naxis = 0;
    long naxes[9] = {1,1,1,1,1,1,1,1,1};
    fits_get_img_param(input, 9, &bitpix, &naxis, naxes, &status);
    if (status) throw std::runtime_error("cannot read logical image parameters: " + fitsDiagnostic(status));
    if (naxis != 2) throw std::runtime_error("initializer only extracts two-dimensional HDUs");
    fitsfile* output = nullptr;
    try {
        const std::string name = "!" + output_path.string();
        fits_create_file(&output, name.c_str(), &status);
        if (status == 0) fits_create_img(output, bitpix, naxis, naxes, &status);
        if (status) throw std::runtime_error("cannot create uncompressed output: " + fitsDiagnostic(status));
        int count = 0;
        fits_get_hdrspace(input, &count, nullptr, &status);
        if (status) {
            const auto msg = fitsDiagnostic(status);
            if (!relaxed) throw HeaderFailure("cannot enumerate header: " + msg);
            throw std::runtime_error("cannot enumerate repaired header: " + msg);
        }
        for (int i = 1; i <= count; ++i) {
            char card[FLEN_CARD] = {};
            int rs = 0;
            fits_read_record(input, i, card, &rs);
            if (rs) {
                const auto msg = fitsDiagnostic(rs);
                if (!relaxed) throw HeaderFailure("cannot read header record: " + msg);
                std::cerr << "[init] " << source << " HDU=" << physical_hdu
                          << " DROP unreadable optional record " << i << ": " << msg << '\n';
                continue;
            }
            const std::string c(card, std::strlen(card));
            const std::string key = trim(c.substr(0, std::min<std::size_t>(8,c.size())));
            if (relaxed) {
                std::string padded = c;
                padded.resize(80, ' ');
                if (invalidOptionalCard(padded)) {
                    if (critical(key)) throw std::runtime_error("malformed critical header card " + key);
                    std::cerr << "[init] " << source << " HDU=" << physical_hdu
                              << " DROP malformed optional header card " << key << '\n';
                    continue;
                }
            }
            const int cls = fits_get_keyclass(card);
            if (cls <= TYP_CMPRS_KEY || cls == TYP_CKSUM_KEY) continue;
            int ws = 0;
            fits_write_record(output, card, &ws);
            if (!ws) continue;
            const auto msg = fitsDiagnostic(ws);
            if (!relaxed) throw HeaderFailure("cannot copy FITS header cards: " + msg);
            if (critical(key)) throw std::runtime_error("cannot copy required header " + key + ": " + msg);
            std::cerr << "[init] " << source << " HDU=" << physical_hdu
                      << " DROP rejected optional header " << key << ": " << msg << '\n';
        }
        const auto [datatype, bytes_per_pixel] = imageDatatype(bitpix);
        LONGLONG total = 1;
        for (int a = 0; a < naxis; ++a) {
            if (naxes[a] <= 0 || total > LLONG_MAX / naxes[a])
                throw std::runtime_error("invalid FITS image axis/pixel count");
            total *= naxes[a];
        }
        const LONGLONG chunk = std::max<LONGLONG>(1, kCopyBufferBytes / bytes_per_pixel);
        const auto bufbytes = static_cast<std::size_t>(chunk) * bytes_per_pixel;
        std::vector<double> buffer((bufbytes + sizeof(double)-1) / sizeof(double));
        double nullvalue = 0;
        int anynull = 0;
        fits_set_bscale(input, 1.0, 0.0, &status);
        if (status == 0) fits_set_bscale(output, 1.0, 0.0, &status);
        for (LONGLONG first = 1; first <= total && status == 0;) {
            const LONGLONG n = std::min(chunk, total - first + 1);
            fits_read_img(input, datatype, first, n, &nullvalue, buffer.data(), &anynull, &status);
            if (status == 0) fits_write_img(output, datatype, first, n, buffer.data(), &status);
            first += n;
        }
        int cs = 0;
        fits_close_file(output, &cs);
        output = nullptr;
        if (status == 0) status = cs;
        if (status) throw std::runtime_error("cannot stream FITS pixels: " + fitsDiagnostic(status));
    } catch (...) { closeFits(output); throw; }
}

// Scan only raw 80-byte cards to avoid calling CFITSIO's vulnerable compressed
// HDU counter on untrusted metadata. No original image pixels are examined.
bool integerFrom(const std::map<std::string,std::string>& h, const std::string& key,
                 long long& value) {
    auto it = h.find(key);
    if (it == h.end() || !validInteger(it->second)) return false;
    try {
        std::size_t n = 0;
        value = std::stoll(it->second, &n);
        return n == it->second.size();
    } catch (...) { return false; }
}
bool checkedMul(std::uint64_t x, std::uint64_t y, std::uint64_t& out) {
    if (y && x > UINT64_MAX/y) return false;
    out = x*y;
    return true;
}
struct ScannedHeader {
    std::map<std::string,std::string> values;
    std::set<std::string> string_values;
    std::uint64_t bytes = 0;
    bool unsafe = false, sanitize = false;
};
ScannedHeader readHeader(std::ifstream& in, std::uint64_t offset, bool primary) {
    ScannedHeader h;
    in.clear();
    in.seekg(static_cast<std::streamoff>(offset));
    if (!in) throw std::runtime_error("cannot seek FITS header");
    bool ended = false;
    std::size_t idx = 0;
    for (unsigned b = 0; b < kMaxHeaderBlocks && !ended; ++b) {
        std::array<char, 2880> block{};
        in.read(block.data(), block.size());
        if (!in) throw std::runtime_error("truncated FITS header");
        for (int j = 0; j < 36; ++j, ++idx) {
            const std::string c(block.data()+j*80,80);
            const std::string key = trim(c.substr(0,8));
            if (idx == 0 && key != (primary ? "SIMPLE" : "XTENSION")) h.unsafe = true;
            if (key == "END") {
                if (c != std::string("END") + std::string(77, ' ')) h.sanitize = true;
                ended = true;
                break;
            }
            if (invalidOptionalCard(c)) {
                if (critical(key) || idx == 0) h.unsafe = true;
                else h.sanitize = true;
            }
            if (key.empty() || c[8] != '=') continue;
            std::string v = trim(c.substr(10));
            if (!v.empty() && v.front() == '\'') {
                h.string_values.insert(key);
                std::size_t close = v.find('\'',1);
                if (close != std::string::npos) v = trim(v.substr(1,close-1));
            } else {
                v = trim(v.substr(0,v.find('/')));
            }
            if (h.values.count(key) && critical(key)) h.unsafe = true;
            if (!h.values.count(key)) h.values[key] = v;
        }
        if (ended) h.bytes = (b+1)*kFitsBlock;
    }
    if (!ended) throw std::runtime_error("missing FITS END card");
    return h;
}
std::vector<RawHdu> scanRaw(const fs::path& source, bool& incomplete) {
    const auto size = fs::file_size(source);
    std::ifstream in(source, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open FITS source");
    std::uint64_t pos = 0;
    std::vector<RawHdu> hdus;
    incomplete = false;
    while (pos < size) {
      try {
        auto head = readHeader(in, pos, pos == 0);
        const auto& h = head.values;
        long long bitpix=0,naxis=0,pcount=0,gcount=1;
        if (!integerFrom(h,"BITPIX",bitpix) || !integerFrom(h,"NAXIS",naxis) ||
            naxis < 0 || naxis > 9 ||
            (bitpix!=8 && bitpix!=16 && bitpix!=32 && bitpix!=64 && bitpix!=-32 && bitpix!=-64))
            throw std::runtime_error("invalid physical BITPIX/NAXIS at byte " + std::to_string(pos));
        if (h.count("PCOUNT") && !integerFrom(h,"PCOUNT",pcount))
            throw std::runtime_error("invalid physical PCOUNT");
        if (h.count("GCOUNT") && !integerFrom(h,"GCOUNT",gcount))
            throw std::runtime_error("invalid physical GCOUNT");
        if (pcount < 0 || gcount <= 0) throw std::runtime_error("invalid physical PCOUNT/GCOUNT");
        std::uint64_t elements = naxis == 0 ? 0 : 1;
        for (int i = 1; i <= naxis; ++i) {
            long long n = 0;
            if (!integerFrom(h,"NAXIS"+std::to_string(i),n) || n < 0 ||
                !checkedMul(elements,static_cast<std::uint64_t>(n),elements))
                throw std::runtime_error("unreadable physical NAXISn at byte " + std::to_string(pos));
        }
        const bool table = h.count("XTENSION") &&
            (h.at("XTENSION") == "BINTABLE" || h.at("XTENSION") == "TABLE");
        std::uint64_t nbytes=0, tmp=0;
        if (elements > UINT64_MAX-static_cast<std::uint64_t>(pcount) ||
            !checkedMul(elements+static_cast<std::uint64_t>(pcount),
                        static_cast<std::uint64_t>(gcount),tmp) ||
            (!table && !checkedMul(tmp,static_cast<std::uint64_t>(std::abs(bitpix)/8),nbytes)))
            throw std::runtime_error("FITS physical byte-count overflow");
        if (table) nbytes=tmp;
        if (nbytes > UINT64_MAX-(kFitsBlock-1)) throw std::runtime_error("FITS padding overflow");
        const auto padded=(nbytes+kFitsBlock-1)/kFitsBlock*kFitsBlock;
        if (pos > size || head.bytes > size-pos || padded > size-pos-head.bytes)
            throw std::runtime_error("truncated FITS extension");
        RawHdu r;
        r.hdu=static_cast<int>(hdus.size())+1;
        r.offset=pos; r.header_length=head.bytes; r.total_length=head.bytes+padded;
        r.unsafe=head.unsafe; r.sanitize=head.sanitize;
        r.compressed=table && h.count("ZIMAGE") && h.at("ZIMAGE")=="T";
        long long dim=naxis;
        if (r.compressed) {
            r.image=true;
            if (!integerFrom(h,"ZNAXIS",dim) || dim!=2) r.unsafe=true;
            if (dim >= 0 && dim <= 9) {
                for (int a=1;a<=dim;++a) {
                    long long length=0, tile=0;
                    const auto k="ZNAXIS"+std::to_string(a);
                    if (!integerFrom(h,k,length) || length <= 0) { r.unsafe=true; continue; }
                    const auto tk="ZTILE"+std::to_string(a);
                    if (h.count(tk) && (!integerFrom(h,tk,tile) || tile<=0)) r.unsafe=true;
                    // Absent ZTILE is a legitimate CFITSIO default.
                }
            }
        } else {
            r.image=!table && naxis==2 &&
                    (pos==0 || (h.count("XTENSION") && h.at("XTENSION")=="IMAGE"));
        }
        r.naxis=static_cast<int>(dim);
        long long ccd=0;
        if (!head.string_values.count(InitConfig::CCDNUM_KEYWORD) &&
            integerFrom(h,InitConfig::CCDNUM_KEYWORD,ccd) &&
            ccd>=INT_MIN && ccd<=INT_MAX) {
            r.has_ccdnum=true; r.ccdnum=static_cast<int>(ccd);
        }
        hdus.push_back(r);
        pos+=r.total_length;
      } catch (const std::exception& e) {
        if (hdus.empty()) throw;
        incomplete=true;
        std::cerr << "[init] " << source << " physical scan ended after " << hdus.size()
                  << " HDUs: " << e.what() << '\n';
        break;
      }
    }
    return hdus;
}

void copyRange(const fs::path& path, std::ofstream& out, std::uint64_t offset,
               std::uint64_t length) {
    std::ifstream in(path,std::ios::binary);
    if (!in) throw std::runtime_error("cannot reopen source file");
    in.seekg(static_cast<std::streamoff>(offset));
    if (!in) throw std::runtime_error("source seek failed");
    std::array<char,1024*1024> buf{};
    while (length) {
        auto n=static_cast<std::streamsize>(std::min<std::uint64_t>(length,buf.size()));
        in.read(buf.data(),n);
        if (in.gcount()!=n) throw std::runtime_error("truncated temporary slice input");
        out.write(buf.data(),n);
        if (!out) throw std::runtime_error("temporary FITS write failed");
        length-=static_cast<std::uint64_t>(n);
    }
}
void sanitizedHdu(const fs::path& source,std::ofstream& out,const PlannedImage& plan) {
    std::ifstream in(source,std::ios::binary);
    if (!in) throw std::runtime_error("cannot read damaged header");
    in.seekg(static_cast<std::streamoff>(plan.offset));
    if (!in) throw std::runtime_error("cannot seek damaged header");
    std::vector<std::array<char,80>> cards;
    bool ended=false;
    std::size_t index=0;
    for (std::uint64_t off=0;off<plan.header_length && !ended;off+=kFitsBlock) {
        std::array<char,2880> block{};
        in.read(block.data(),block.size());
        if (!in) throw std::runtime_error("short source header");
        for (int j=0;j<36;++j,++index) {
            std::string card(block.data()+j*80,80);
            const auto key=trim(card.substr(0,8));
            if (index==0 && key!=(plan.hdu==1?"SIMPLE":"XTENSION"))
                throw std::runtime_error("damaged first mandatory HDU header card");
            if (key=="END") { ended=true; break; }
            if (invalidOptionalCard(card)) {
                if (critical(key) || index==0)
                    throw std::runtime_error("invalid critical FITS card "+key);
                std::cerr << "[init] " << source << " HDU=" << plan.hdu
                          << " DROP malformed optional card " << index+1
                          << " key=" << key << '\n';
                continue;
            }
            std::array<char,80> c{};
            std::copy(card.begin(),card.end(),c.begin());
            cards.push_back(c);
        }
    }
    if (!ended) throw std::runtime_error("temporary FITS missing END");
    std::array<char,80> last{};
    last.fill(' ');
    std::memcpy(last.data(),"END",3);
    cards.push_back(last);
    for (std::size_t base=0;base<cards.size();base+=36) {
        std::array<char,2880> block{};
        block.fill(' ');
        for (std::size_t j=0;j<36 && base+j<cards.size();++j)
            std::copy(cards[base+j].begin(),cards[base+j].end(),block.data()+j*80);
        out.write(block.data(),block.size());
        if (!out) throw std::runtime_error("sanitized FITS header write failed");
    }
    copyRange(source,out,plan.offset+plan.header_length,plan.total_length-plan.header_length);
}
void emptyPrimary(std::ofstream& out) {
    std::array<char,2880> block{};
    block.fill(' ');
    const char* cards[]={"SIMPLE  =                    T","BITPIX  =                    8",
                         "NAXIS   =                    0","EXTEND  =                    T","END"};
    for (std::size_t i=0;i<5;++i) std::memcpy(block.data()+i*80,cards[i],std::strlen(cards[i]));
    out.write(block.data(),block.size());
}

std::vector<PlannedImage> makeRawPlans(const fs::path& source,ProductKind kind,
         const fs::path& target,const fs::path& stage,
         const std::vector<RawHdu>& hdus,int& skipped) {
    const std::string stem=kind==ProductKind::Science?archiveStem(source):dqOutputStem(source);
    std::set<int> numbers;
    std::vector<PlannedImage> plans;
    for (const auto& h:hdus) {
        if (!h.image) continue;
        if (h.unsafe || h.naxis!=2 || !h.has_ccdnum) {
            ++skipped;
            std::cerr << "[init] " << source << " HDU=" << h.hdu
                      << " SKIP unsafe structure, axes or native-integer CCDNUM\n";
            continue;
        }
        if (!numbers.insert(h.ccdnum).second) {
            ++skipped;
            std::cerr << "[init] " << source << " HDU=" << h.hdu
                      << " SKIP duplicate CCDNUM=" << h.ccdnum << '\n';
            continue;
        }
        auto filename=stem+"_"+std::to_string(h.ccdnum)+".fits";
        plans.push_back({h.hdu,h.ccdnum,target/filename,stage/filename,true,h.sanitize,
                         h.offset,h.header_length,h.total_length});
    }
    return plans;
}
std::vector<PlannedImage> makeNormalPlans(fitsfile* input,const fs::path& source,
       ProductKind kind,const fs::path& target,const fs::path& stage,int& skipped) {
    int count=0,status=0;
    fits_get_num_hdus(input,&count,&status);
    if (status) throw std::runtime_error("cannot count FITS HDUs: "+fitsDiagnostic(status));
    const std::string stem=kind==ProductKind::Science?archiveStem(source):dqOutputStem(source);
    std::set<int> seen;
    std::vector<PlannedImage> plans;
    for (int h=1;h<=count;++h) {
        int type=0,s=0,naxis=0;
        fits_movabs_hdu(input,h,&type,&s);
        if (s || type!=IMAGE_HDU) { fits_clear_errmsg(); continue; }
        fits_get_img_dim(input,&naxis,&s);
        if (s || naxis!=2) { fits_clear_errmsg(); continue; }
        if (!nativeKey(input,InitConfig::CCDNUM_KEYWORD,true)) {
            ++skipped;
            std::cerr << "[init] " << source << " HDU=" << h
                      << " SKIP CCDNUM is not native FITS integer\n";
            continue;
        }
        int ccd=0;
        fits_read_key(input,TINT,InitConfig::CCDNUM_KEYWORD,&ccd,nullptr,&s);
        if (s) {
            ++skipped;
            std::cerr << "[init] " << source << " HDU=" << h
                      << " SKIP CCDNUM conversion: " << fitsDiagnostic(s) << '\n';
            continue;
        }
        if (!seen.insert(ccd).second) {
            ++skipped;
            std::cerr << "[init] " << source << " HDU=" << h
                      << " SKIP duplicate CCDNUM=" << ccd << '\n';
            continue;
        }
        auto filename=stem+"_"+std::to_string(ccd)+".fits";
        plans.push_back({h,ccd,target/filename,stage/filename});
    }
    return plans;
}

bool validateOutput(const PlannedImage& plan) {
    fitsfile* f=nullptr;
    int status=0,naxis=0;
    fits_open_file(&f,plan.final_path.string().c_str(),READONLY,&status);
    if (status) { fitsDiagnostic(status); closeFits(f); return false; }
    fits_get_img_dim(f,&naxis,&status);
    int compressed=status==0?fits_is_compressed_image(f,&status):0;
    if (status || naxis!=2 || compressed || !nativeKey(f,InitConfig::CCDNUM_KEYWORD,true)) {
        closeFits(f); return false;
    }
    int ccd=0;
    fits_read_key(f,TINT,InitConfig::CCDNUM_KEYWORD,&ccd,nullptr,&status);
    closeFits(f);
    return status==0 && ccd==plan.ccdnum;
}
void extractOne(const fs::path& source,const PlannedImage& plan,ProductKind kind,
                const fs::path& stage,fitsfile* original) {
    fitsfile* f=original;
    fs::path slice;
    if (plan.isolated) {
        slice=stage/(".fq_raw_hdu_"+std::to_string(plan.hdu)+".fits");
        {
            std::ofstream out(slice,std::ios::binary|std::ios::trunc);
            if (!out) throw std::runtime_error("cannot create temporary isolated HDU");
            if (plan.hdu!=1) emptyPrimary(out);
            if (plan.sanitize) sanitizedHdu(source,out,plan);
            else copyRange(source,out,plan.offset,plan.total_length);
            out.close();
            if (!out) throw std::runtime_error("cannot finalize temporary isolated HDU");
        }
        int status=0,type=0;
        fits_open_file(&f,slice.string().c_str(),READONLY,&status);
        if (!status) fits_movabs_hdu(f,plan.hdu==1?1:2,&type,&status);
        if (status || type!=IMAGE_HDU) {
            auto msg=fitsDiagnostic(status);
            closeFits(f);
            std::error_code ec; fs::remove(slice,ec);
            throw std::runtime_error("cannot open isolated HDU: "+msg);
        }
    } else {
        int type=0,s=0;
        fits_movabs_hdu(f,plan.hdu,&type,&s);
        if (s || type!=IMAGE_HDU)
            throw std::runtime_error("cannot revisit planned image HDU: "+fitsDiagnostic(s));
    }
    try {
        validateRequired(f,kind);
        try { copyCurrentImage(f,plan.staged_path,plan.sanitize,source,plan.hdu); }
        catch (const HeaderFailure& e) {
            std::cerr << "[init] " << source << " HDU=" << plan.hdu
                      << " retry optional-header copy: " << e.what() << '\n';
            fits_clear_errmsg();
            copyCurrentImage(f,plan.staged_path,true,source,plan.hdu);
        }
        if (plan.isolated) {
            closeFits(f);
            std::error_code ec; fs::remove(slice,ec);
        }
    } catch (...) {
        if (plan.isolated) {
            closeFits(f);
            std::error_code ec; fs::remove(slice,ec);
        }
        throw;
    }
}
} // namespace

std::string archiveStem(const fs::path& source) {
    const auto name=source.filename().string();
    const auto n=std::strlen(InitConfig::ARCHIVE_SUFFIX);
    if (name.size()<=n || name.compare(name.size()-n,n,InitConfig::ARCHIVE_SUFFIX)!=0)
        throw std::runtime_error("archive suffix mismatch: "+source.string());
    return name.substr(0,name.size()-n);
}
std::string dqOutputStem(const fs::path& source) {
    auto stem=archiveStem(source);
    replaceAll(stem,InitConfig::DQ_STEM_REPLACE_FROM,InitConfig::DQ_STEM_REPLACE_TO);
    return stem;
}
ExtractionResult extractArchive(const fs::path& source,ProductKind kind,
      const fs::path& final_directory,const fs::path& staging_directory,
      ExistingPolicy policy) {
    ExtractionResult result;
    fitsfile* input=nullptr;
    try {
        bool incomplete=false, salvage=false;
        const auto hdus=scanRaw(source,incomplete);
        salvage=incomplete;
        for (const auto& h:hdus) salvage|=h.unsafe || h.sanitize;
        int skipped=incomplete?1:0;
        std::vector<PlannedImage> plans;
        if (!salvage) {
            int s=0;
            fits_open_file(&input,source.string().c_str(),READONLY,&s);
            if (s) throw std::runtime_error("cannot open archive: "+fitsDiagnostic(s));
            try { plans=makeNormalPlans(input,source,kind,final_directory,staging_directory,skipped); }
            catch (const std::exception& e) {
                closeFits(input);
                skipped=incomplete?1:0;
                salvage=true;
                std::cerr << "[init] " << source << " fallback after HDU discovery error: "
                          << e.what() << '\n';
            }
        }
        if (salvage) {
            closeFits(input);
            std::cerr << "[init] " << source
                      << " malformed or dangerous FITS header detected; salvaging other HDUs\n";
            plans=makeRawPlans(source,kind,final_directory,staging_directory,hdus,skipped);
        }
        if (plans.empty()) throw std::runtime_error("archive contains no extractable two-dimensional image HDUs");
        bool existing=false;
        for (const auto& p:plans) existing|=fs::exists(p.final_path);
        if (existing && policy==ExistingPolicy::Fail)
            throw std::runtime_error("output already exists: "+final_directory.string());
        std::vector<PlannedImage> pending;
        for (const auto& p:plans) {
            if (policy==ExistingPolicy::Resume && fs::exists(p.final_path) && validateOutput(p)) {
                result.output_paths.push_back(p.final_path);
            } else pending.push_back(p);
        }
        if (pending.empty()) {
            result.success=!result.output_paths.empty();
            result.resumed=result.success;
            result.skipped_hdus=skipped;
            closeFits(input);
            return result;
        }
        fs::create_directories(staging_directory);
        fs::create_directories(final_directory);
        std::vector<PlannedImage> successful;
        for (const auto& p:pending) {
            try { extractOne(source,p,kind,staging_directory,input); successful.push_back(p); }
            catch (const std::exception& e) {
                ++skipped;
                std::cerr << "[init] " << source << " HDU=" << p.hdu
                          << " SKIP: " << e.what() << '\n';
                fits_clear_errmsg();
                std::error_code ec; fs::remove(p.staged_path,ec);
            }
        }
        closeFits(input);
        for (const auto& p:successful) {
            if (std::rename(p.staged_path.c_str(),p.final_path.c_str()) != 0) {
                ++skipped;
                std::cerr << "[init] " << source << " HDU=" << p.hdu
                          << " SKIP commit " << p.final_path << ": " << std::strerror(errno) << '\n';
                continue;
            }
            result.output_paths.push_back(p.final_path);
        }
        result.skipped_hdus=skipped;
        result.success=!result.output_paths.empty();
        if (!result.success) result.error="all planned HDUs failed extraction";
        std::error_code ec; fs::remove_all(staging_directory,ec);
        return result;
    } catch (const std::exception& e) {
        closeFits(input);
        result.error=e.what();
        std::error_code ec; fs::remove_all(staging_directory,ec);
        return result;
    }
}
} // namespace fqinit
