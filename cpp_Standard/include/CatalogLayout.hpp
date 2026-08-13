#ifndef CATALOG_LAYOUT_HPP
#define CATALOG_LAYOUT_HPP

#include "RuntimeConfig.hpp"

#include <cstddef>
#include <optional>
#include <string>

namespace PipelineCatalog {

// ==========================================
// Structure: Describe named external-catalog field positions
// Method: Keep RA/Dec/ZP mandatory while representing independently optional
//         magnitude bands as effective zero-based positions when available.
// ==========================================
struct ExternalFields {
    std::size_t ra = 0;
    std::size_t dec = 0;
    std::optional<std::size_t> mag_g;
    std::optional<std::size_t> mag_r;
    std::optional<std::size_t> mag_i;
    std::optional<std::size_t> mag_z;
    std::optional<std::size_t> mag_y;
    std::size_t zp = 0;
};

// ==========================================
// Structure: Describe process_main fields in complete _all.cat rows
// Method: Store zero-based absolute positions after the runtime external prefix and CCD field.
// ==========================================
struct SourceFields {
    std::size_t polychi2 = 0;
    std::size_t pixx = 0;
    std::size_t pixy = 0;
    std::size_t sig = 0;
    std::size_t star = 0;
    std::size_t peak = 0;
    std::size_t imax = 0;
    std::size_t jmax = 0;
    std::size_t h_flux = 0;
    std::size_t h_area = 0;
    std::size_t flag = 0;
    std::size_t psf = 0;
    std::size_t snr_f = 0;
    std::size_t ra = 0;
    std::size_t dec = 0;
    std::size_t gf1 = 0;
    std::size_t gf2 = 0;
    std::size_t g1 = 0;
    std::size_t g2 = 0;
    std::size_t de = 0;
    std::size_t h1 = 0;
    std::size_t h2 = 0;
    std::size_t cos2 = 0;
    std::size_t sin2 = 0;
    std::size_t parity = 0;
    std::size_t galsizeT = 0;
    std::size_t psfsizeT = 0;
    std::size_t delta_chi2 = 0;
    std::size_t orth_ext = 0;
    std::size_t chi2 = 0;
};

// ==========================================
// Structure: Define the shared runtime schema for every downstream catalog phase
// Method: Combine effective external fields, CCD position, process_main suffix, and row widths.
// ==========================================
struct CatalogLayout {
    std::size_t external_columns = 0;
    ExternalFields external;
    std::size_t ccd = 0;
    std::size_t source_base = 0;
    std::size_t source_columns = 0;
    SourceFields source;
    std::size_t all_columns = 0;
};

// ==========================================
// Function: Resolve the pipeline catalog schema once from runtime options
// Method: Map mandatory fields and available optional magnitudes through the
//         projection, then derive every downstream row offset.
// ==========================================
bool resolveCatalogLayout(const RuntimeConfig& config,
                          CatalogLayout& layout,
                          std::string& error);

// ==========================================
// Structure: RearrCatalogSchema
// Method: Describe only the complete-row width and sky columns needed by
//         process_rearr for either external or internal Stage-9 output.
// ==========================================
struct RearrCatalogSchema {
    std::size_t all_columns = 0;
    std::size_t ra_column = 0;
    std::size_t dec_column = 0;
};

// ==========================================
// Function: resolveRearrCatalogSchema
// Method: Derive rearr fields from the external layout when ext_cat=1 or from
//         the fixed internal `[CCD_NUM]+source` schema when ext_cat=0.
// ==========================================
bool resolveRearrCatalogSchema(const RuntimeConfig& config,
                               const CatalogLayout* external_layout,
                               RearrCatalogSchema& schema,
                               std::string& error);

// ==========================================
// Function: Describe one resolved catalog layout
// Method: Render the effective widths and key zero-based offsets for startup diagnostics.
// ==========================================
std::string describeCatalogLayout(const CatalogLayout& layout);

}  // namespace PipelineCatalog

#endif  // CATALOG_LAYOUT_HPP
