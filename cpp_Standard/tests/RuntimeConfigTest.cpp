#include "RuntimeConfig.hpp"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
int assertions=0, failures=0;
void expect(bool ok,const std::string& message) {
    ++assertions;
    if (!ok) { ++failures; std::cerr<<"FAIL: "<<message<<'\n'; }
}
template<typename Lens>
auto checkAlias(const RuntimeConfig& c,const Lens& lens,int) -> decltype((void)lens.source_cat) {
    expect(c.extcat.output_directory=="/catalog/a" && lens.source_cat=="/catalog/b",
           "Lite independent producer/consumer");
}
template<typename Lens>
void checkAlias(const RuntimeConfig& c,const Lens&,long) {
    expect(c.extcat.output_directory=="/catalog/b","Standard source_cat alias");
}
bool scan(std::vector<std::string> words,std::string& path) {
    std::vector<char*> argv;
    for (auto& word:words) argv.push_back(&word[0]);
    bool help=true; std::string error;
    bool ok=findRuntimeConfigPath(static_cast<int>(argv.size()),argv.data(),path,help,error);
    expect(!help,"no help CLI flag");
    return ok;
}
}
int main(int argc,char* argv[]) {
    std::string path,error;
    setenv("FQ_PIPELINE_CONFIG","/must/not/be/read.ini",1);
    expect(scan({"pipeline"},path) && path.empty(),"no arguments retain compiled configuration");
    expect(scan({"pipeline","--config","/data/my file.ini"},path) && path=="/data/my file.ini","explicit config separated");
    expect(scan({"pipeline","--config=/data/file.ini"},path) && path=="/data/file.ini","explicit config equals");
    for (const auto& args:std::vector<std::vector<std::string>>{
         {"pipeline","--config"},{"pipeline","--config="},{"pipeline","--config",""},
         {"pipeline","--config","a","--config","b"}, {"pipeline","a.list"},
         {"pipeline","--config","a","trailing"}})
        expect(!scan(args,path),"invalid argument sequence is rejected");
    for (const std::string option:{"--help","-h","--run-astrocat","--astrocat-input","--astrocat-output",
        "--astrocat-add-header","--astrocat-existing","--run-extcat","--run-init","--run-main",
        "--run-rearr","--run-fd","--extcat-input","--extcat-output","--extcat-contains",
        "--extcat-recursive","--extcat-delimiter","--extcat-header","--extcat-columns",
        "--extcat-ra-column","--extcat-dec-column","--extcat-zp-column","--extcat-chunk-mib",
        "--extcat-malformed","--extcat-existing","--science-root","--dq-root","--output-root",
        "--dataset","--target","--prefix","--contains","--existing","--f77-max-path","--expo-list",
        "--rearr-output-dir","--rearr-output-base","--rearr-list-name","--rearr-list-dir",
        "--fd-expo-list","--fd-output-dir","--fd-output-base"})
        expect(!scan({"pipeline",option,"value"},path),"removed CLI rejected: "+option);
    auto defaults=makeDefaultRuntimeConfig();
    auto c=defaults;
    expect(applyRuntimeConfigText("[process]\n# run_process_init =\n[lensing]\n# chipnx =\n", "comments", c,error),error);
    expect(c.process.run_process_init==defaults.process.run_process_init && c.lensing.chipnx==defaults.lensing.chipnx,
           "commented blanks inherit defaults");
    expect(applyRuntimeConfigText("[init]\nscience_root = /new/path\n", "paths",c,error),error);
    expect(c.init.science_root=="/new/path" && c.lensing.chipnx==defaults.lensing.chipnx,"paths do not replace scientific defaults");
    expect(applyRuntimeConfigText("[init]\nscience_root =\n", "empty",c,error),error);
    expect(c.init.science_root.empty(),"active empty string clears, not inherits");
    c=defaults;
    expect(!applyRuntimeConfigText("[init]\nscience_root=/changed\n[process]\nrun_process_init=\n","bad",c,error),"empty bool fails");
    expect(c.init.science_root==defaults.init.science_root,"INI application is transactional");
    for (const std::string entry:{"chipnx=", "chipnx=no", "pixel_size=", "pixel_size=nan", "pixel_size=1abc"})
        expect(!applyRuntimeConfigText("[lensing]\n"+entry+"\n","bad scalar",c,error),"malformed scalar rejected");
    expect(applyRuntimeConfigText("[init]\ncontains=\n","empty list",c,error) && c.init.contains.empty(),"empty tokens clear list");
    expect(!applyRuntimeConfigText("[extcat]\ninput_columns=\n","empty projection",c,error),"projection requires values");
    expect(applyRuntimeConfigText("[process]\nexpo_list=\"/data/a # b.list\" ; comment\n", "quoted",c,error),error);
    expect(c.process.expo_list=="/data/a # b.list" && c.external_expo_list_supplied,"INI exposure list replaces old explicit CLI semantics");
    expect(applyRuntimeConfigText("[extcat]\nuse_explicit_columns=true\ninput_columns=17,5,6\nuse_explicit_coordinate_columns=true\nra_column=5\ndec_column=6\n", "columns",c,error),error);
    expect(c.extcat.use_explicit_columns && c.extcat.use_explicit_coordinate_columns
        && c.extcat.input_columns_one_based==std::vector<std::size_t>({17,5,6}),"projection flags preserved in INI");
    expect(applyRuntimeConfigText("[init]\ndatasets=g:c4d_,r:other_\ncontains=v1,v2\n", "datasets",c,error),error);
    expect(c.init.datasets.size()==2 && c.init.datasets[1].target=="r" && c.init.datasets[1].prefix=="other_",
        "dataset list replaces target/prefix/repeatable dataset CLI");
    expect(!applyRuntimeConfigText("[init]\ndatasets=bad\n","dataset error",c,error),"invalid dataset rejected");
    expect(!applyRuntimeConfigText("[unknown]\nx=1\n","section",c,error),"unknown section rejected");
    expect(!applyRuntimeConfigText("[init]\nunknown=1\n","key",c,error),"unknown key rejected");
    c=defaults;
    expect(applyRuntimeConfigText("[extcat]\noutput_directory=/catalog/a\n[lensing]\nsource_cat=/catalog/b\n", "alias",c,error),error);
    checkAlias(c,c.lensing,0);
    c=defaults;
    expect(applyRuntimeConfigText("[lensing]\nccd_split=99\n","domain",c,error),error);
    expect(!validateRuntimeConfig(c,error),"invalid branch domain fails validation");
    if (argc==2) {
        std::ifstream f(argv[1]); std::ostringstream buffer; buffer<<f.rdbuf();
        expect(f.is_open(),"generated template readable");
        c=defaults;
        expect(applyRuntimeConfigText(buffer.str(),argv[1],c,error),"generated INI accepted by production parser: "+error);
        expect(c.process.run_process_init==defaults.process.run_process_init
            && c.process.run_process_main==defaults.process.run_process_main
            && c.process.run_process_fd==defaults.process.run_process_fd
            && c.lensing.process_stage==defaults.lensing.process_stage
            && c.lensing.chipnx==defaults.lensing.chipnx,"generated template preserves compiled science settings");
    }
    c=defaults;
    expect(RuntimeConfigStore::initialize(c,error),"freeze runtime configuration");
    expect(RuntimeConfigStore::isInitialized(),"store initialized");
    expect(!RuntimeConfigStore::initialize(c,error),"store cannot be replaced");
    std::cout<<assertions<<" assertions, "<<failures<<" failures\n";
    return failures?1:0;
}
