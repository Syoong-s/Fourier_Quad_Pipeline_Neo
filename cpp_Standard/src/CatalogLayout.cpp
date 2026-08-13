#include "CatalogLayout.hpp"

#include "LensingConfig.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace PipelineCatalog {
namespace {

// ==========================================
// Function: Validate an explicit external-catalog projection
// Method: Require positive, unique raw one-based indices so every field has one output position.
// ==========================================
bool validateProjection(const RuntimeConfig& config,
                        std::string& error) {
    const ExtCatRuntimeConfig& extcat = config.extcat;
    if (!extcat.use_explicit_columns) {
        return true;
    }
    if (extcat.input_columns_one_based.empty()) {
        error = "external-catalog explicit projection must not be empty";
        return false;
    }

    std::set<std::size_t> unique_columns;
    for (const std::size_t raw_column :
         extcat.input_columns_one_based) {
        if (raw_column == 0) {
            error = "external-catalog explicit projection columns must be positive one-based indices";
            return false;
        }
        if (!unique_columns.insert(raw_column).second) {
            error = "external-catalog explicit projection repeats raw column "
                    + std::to_string(raw_column);
            return false;
        }
    }
    return true;
}

// ==========================================
// Function: Resolve one mandatory raw field to its effective output position
// Method: Require a positive raw identity and preserve pass-through order or
//         locate that identity in the explicit projection.
// ==========================================
bool resolveRequiredProjectedColumn(
    std::size_t raw_column_one_based,
    const std::string& field_name,
    const RuntimeConfig& config,
    std::size_t& effective_zero_based,
    std::string& error) {
    if (raw_column_one_based == 0) {
        error = "external-catalog " + field_name
                + " raw column must be a positive one-based index";
        return false;
    }
    const ExtCatRuntimeConfig& extcat = config.extcat;
    if (!extcat.use_explicit_columns) {
        effective_zero_based = raw_column_one_based - 1;
        return true;
    }

    const std::vector<std::size_t>& projection =
        extcat.input_columns_one_based;
    const auto match = std::find(projection.begin(), projection.end(),
                                 raw_column_one_based);
    if (match == projection.end()) {
        error = "external-catalog explicit projection omits required "
                + field_name + " raw column "
                + std::to_string(raw_column_one_based);
        return false;
    }
    effective_zero_based = static_cast<std::size_t>(
        std::distance(projection.begin(), match));
    return true;
}

// ==========================================
// Function: Resolve one optional magnitude field
// Method: Interpret raw zero or projection omission as absence; otherwise map
//         the positive raw identity to its effective zero-based position.
// ==========================================
void resolveOptionalProjectedColumn(
    std::size_t raw_column_one_based,
    const RuntimeConfig& config,
    std::optional<std::size_t>& effective_zero_based) {
    effective_zero_based.reset();
    if (raw_column_one_based == 0) {
        return;
    }
    const ExtCatRuntimeConfig& extcat = config.extcat;
    if (!extcat.use_explicit_columns) {
        effective_zero_based = raw_column_one_based - 1;
        return;
    }

    const std::vector<std::size_t>& projection =
        extcat.input_columns_one_based;
    const auto match = std::find(projection.begin(), projection.end(),
                                 raw_column_one_based);
    if (match != projection.end()) {
        effective_zero_based = static_cast<std::size_t>(
            std::distance(projection.begin(), match));
    }
}

// ==========================================
// Function: Validate resolved external field positions
// Method: Require RA/Dec/ZP and every present magnitude to be in range and
//         semantically distinct while permitting absent magnitude bands.
// ==========================================
bool validateExternalFields(const CatalogLayout& layout,
                            std::string& error) {
    const std::array<std::pair<const char*, std::size_t>, 3> required_fields = {{
        {"RA", layout.external.ra},
        {"Dec", layout.external.dec},
        {"ZP", layout.external.zp},
    }};
    const std::array<
        std::pair<const char*, const std::optional<std::size_t>*>, 5>
        optional_fields = {{
            {"mag_g", &layout.external.mag_g},
            {"mag_r", &layout.external.mag_r},
            {"mag_i", &layout.external.mag_i},
            {"mag_z", &layout.external.mag_z},
            {"mag_y", &layout.external.mag_y},
        }};

    std::set<std::size_t> unique_positions;
    const auto validate_position = [&](const char* name,
                                       std::size_t position) -> bool {
        if (position >= layout.external_columns) {
            error = "external-catalog " + std::string(name)
                    + " column lies outside effective width "
                    + std::to_string(layout.external_columns);
            return false;
        }
        if (!unique_positions.insert(position).second) {
            error = "external-catalog field " + std::string(name)
                    + " overlaps another named field";
            return false;
        }
        return true;
    };

    for (const auto& field : required_fields) {
        if (!validate_position(field.first, field.second)) {
            return false;
        }
    }
    for (const auto& field : optional_fields) {
        if (field.second->has_value()
            && !validate_position(field.first, field.second->value())) {
            return false;
        }
    }
    return true;
}

// ==========================================
// Function: Format one optional external-catalog field position
// Method: Print the zero-based effective index when present or "none" when absent.
// ==========================================
std::string describeOptionalColumn(
    const std::optional<std::size_t>& column) {
    return column.has_value() ? std::to_string(column.value()) : "none";
}

}  // namespace

// ==========================================
// Function: Resolve the pipeline catalog schema once from runtime options
// Method: Map mandatory fields and available optional magnitudes through the
//         projection, then derive every downstream row offset.
// ==========================================
bool resolveCatalogLayout(const RuntimeConfig& config,
                          CatalogLayout& layout,
                          std::string& error) {
    if (!validateProjection(config, error)) {
        return false;
    }

    const ExtCatRuntimeConfig& extcat = config.extcat;
    CatalogLayout resolved;
    resolved.external_columns = extcat.use_explicit_columns
        ? extcat.input_columns_one_based.size()
        : extcat.total_columns;
    if (resolved.external_columns == 0) {
        error = "external-catalog effective width must be positive";
        return false;
    }

    if (!resolveRequiredProjectedColumn(
            extcat.ra_column_one_based, "RA", config,
            resolved.external.ra, error)
        || !resolveRequiredProjectedColumn(
            extcat.dec_column_one_based, "Dec", config,
            resolved.external.dec, error)
        || !resolveRequiredProjectedColumn(
            extcat.zp_column_one_based, "ZP", config,
            resolved.external.zp, error)) {
        return false;
    }
    resolveOptionalProjectedColumn(extcat.mag_g_column_one_based,
                                   config, resolved.external.mag_g);
    resolveOptionalProjectedColumn(extcat.mag_r_column_one_based,
                                   config, resolved.external.mag_r);
    resolveOptionalProjectedColumn(extcat.mag_i_column_one_based,
                                   config, resolved.external.mag_i);
    resolveOptionalProjectedColumn(extcat.mag_z_column_one_based,
                                   config, resolved.external.mag_z);
    resolveOptionalProjectedColumn(extcat.mag_y_column_one_based,
                                   config, resolved.external.mag_y);
    if (!validateExternalFields(resolved, error)) {
        return false;
    }

    resolved.source_columns =
        static_cast<std::size_t>(LensingConfig::expo_cat_ncols);
    if (resolved.source_columns == 0
        || resolved.external_columns
               > std::numeric_limits<std::size_t>::max()
                     - 1 - resolved.source_columns) {
        error = "complete _all.cat column count overflows size_t";
        return false;
    }
    resolved.ccd = resolved.external_columns;
    resolved.source_base = resolved.ccd + 1;
    resolved.all_columns = resolved.source_base + resolved.source_columns;

    resolved.source.polychi2 = resolved.source_base + LensingConfig::iid;
    resolved.source.pixx = resolved.source_base + LensingConfig::ipixx;
    resolved.source.pixy = resolved.source_base + LensingConfig::ipixy;
    resolved.source.sig = resolved.source_base + LensingConfig::isig;
    resolved.source.star = resolved.source_base + LensingConfig::istar;
    resolved.source.peak = resolved.source_base + LensingConfig::ipeak;
    resolved.source.imax = resolved.source_base + LensingConfig::i_imax;
    resolved.source.jmax = resolved.source_base + LensingConfig::i_jmax;
    resolved.source.h_flux = resolved.source_base + LensingConfig::ih_flux;
    resolved.source.h_area = resolved.source_base + LensingConfig::ih_area;
    resolved.source.flag = resolved.source_base + LensingConfig::iflag;
    resolved.source.psf = resolved.source_base + LensingConfig::iPSF;
    resolved.source.snr_f = resolved.source_base + LensingConfig::iSNR_F;
    resolved.source.ra = resolved.source_base + LensingConfig::ira;
    resolved.source.dec = resolved.source_base + LensingConfig::idec;
    resolved.source.gf1 = resolved.source_base + LensingConfig::igf1;
    resolved.source.gf2 = resolved.source_base + LensingConfig::igf2;
    resolved.source.g1 = resolved.source_base + LensingConfig::ig1;
    resolved.source.g2 = resolved.source_base + LensingConfig::ig2;
    resolved.source.de = resolved.source_base + LensingConfig::ide;
    resolved.source.h1 = resolved.source_base + LensingConfig::ih1;
    resolved.source.h2 = resolved.source_base + LensingConfig::ih2;
    resolved.source.cos2 = resolved.source_base + LensingConfig::icos2;
    resolved.source.sin2 = resolved.source_base + LensingConfig::isin2;
    resolved.source.parity = resolved.source_base + LensingConfig::iparity;
    resolved.source.galsizeT = resolved.source_base + LensingConfig::igalsizeT;
    resolved.source.psfsizeT = resolved.source_base + LensingConfig::ipsfsizeT;
    resolved.source.delta_chi2 =
        resolved.source_base + LensingConfig::idelta_chi2;
    resolved.source.orth_ext =
        resolved.source_base + LensingConfig::iorth_ext;
    resolved.source.chi2 = resolved.source_base + LensingConfig::ichi2;

    if (resolved.ccd != resolved.external_columns
        || resolved.source_base != resolved.ccd + 1
        || resolved.source.chi2 + 1 != resolved.all_columns) {
        error = "resolved catalog layout is internally inconsistent";
        return false;
    }

    layout = resolved;
    error.clear();
    return true;
}

// ==========================================
// Function: Resolve process_rearr's minimal runtime schema
// Method: Reuse the external layout when enabled or derive the 30-field internal
//         `[CCD_NUM]+28 Stage-7 fields+Chi2` contract directly from fixed indices.
// ==========================================
bool resolveRearrCatalogSchema(const RuntimeConfig& config,
                               const CatalogLayout* external_layout,
                               RearrCatalogSchema& schema,
                               std::string& error) {
    RearrCatalogSchema resolved;
    if (config.lensing.ext_cat == 1) {
        if (external_layout == nullptr) {
            error = "external rearr schema requires a resolved CatalogLayout";
            return false;
        }
        resolved.all_columns = external_layout->all_columns;
        resolved.ra_column = external_layout->external.ra;
        resolved.dec_column = external_layout->external.dec;
    } else {
        resolved.all_columns = 1U
            + static_cast<std::size_t>(LensingConfig::expo_cat_ncols);
        resolved.ra_column = 1U + static_cast<std::size_t>(LensingConfig::ira);
        resolved.dec_column = 1U + static_cast<std::size_t>(LensingConfig::idec);
    }
    if (resolved.all_columns == 0
        || resolved.ra_column >= resolved.all_columns
        || resolved.dec_column >= resolved.all_columns
        || resolved.ra_column == resolved.dec_column) {
        error = "resolved process_rearr catalog schema is inconsistent";
        return false;
    }
    schema = resolved;
    error.clear();
    return true;
}

// ==========================================
// Function: Describe one resolved catalog layout
// Method: Render the effective widths and key zero-based offsets for startup diagnostics.
// ==========================================
std::string describeCatalogLayout(const CatalogLayout& layout) {
    std::ostringstream output;
    output << "Catalog layout: external_columns=" << layout.external_columns
           << " external_ra=" << layout.external.ra
           << " external_dec=" << layout.external.dec
           << " external_mag_g="
           << describeOptionalColumn(layout.external.mag_g)
           << " external_mag_r="
           << describeOptionalColumn(layout.external.mag_r)
           << " external_mag_i="
           << describeOptionalColumn(layout.external.mag_i)
           << " external_mag_z="
           << describeOptionalColumn(layout.external.mag_z)
           << " external_mag_y="
           << describeOptionalColumn(layout.external.mag_y)
           << " external_zp=" << layout.external.zp
           << " ccd=" << layout.ccd
           << " source_base=" << layout.source_base
           << " source_columns=" << layout.source_columns
           << " all_columns=" << layout.all_columns
           << " chi2=" << layout.source.chi2
           << " (indices are zero-based)";
    return output.str();
}

}  // namespace PipelineCatalog
