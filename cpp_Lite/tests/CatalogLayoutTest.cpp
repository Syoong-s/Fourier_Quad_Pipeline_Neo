#include "CatalogLayout.hpp"
#include "process_fd/FDMagnitudeSelector.hpp"
#include "process_fd/ShearCatalogReader.hpp"
#include "process_main/ExternalCatalogReader.hpp"

#include <cstdio>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

// ==========================================
// Function: Record one test expectation
// Method: Count failures and print only violated invariants for concise test output.
// ==========================================
void expect(bool condition, const std::string& message, int& failures) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

// ==========================================
// Function: Create explicit projection options
// Method: Start from production defaults and replace only projection controls.
// ==========================================
ProcessConfig::RuntimeOptions projectionOptions(
    std::vector<std::size_t> columns) {
    ProcessConfig::RuntimeOptions options;
    options.extcat_use_explicit_columns = true;
    options.extcat_input_columns_one_based = std::move(columns);
    return options;
}

// ==========================================
// Function: Build one whitespace-delimited numeric row
// Method: Emit a deterministic sequence with exactly the requested number of columns.
// ==========================================
std::string numericRow(std::size_t columns) {
    std::ostringstream output;
    for (std::size_t column = 0; column < columns; ++column) {
        if (column != 0) {
            output << ' ';
        }
        output << column + 1;
    }
    return output.str();
}

// ==========================================
// Function: Build one numeric row from explicit values
// Method: Preserve the supplied order for synthetic FD reader coverage.
// ==========================================
std::string numericRow(const std::vector<float>& values) {
    std::ostringstream output;
    for (std::size_t column = 0; column < values.size(); ++column) {
        if (column != 0) {
            output << ' ';
        }
        output << values[column];
    }
    return output.str();
}

// ==========================================
// Function: Mark every configured magnitude as absent
// Method: Apply the raw one-based zero sentinel to all five RuntimeOptions fields.
// ==========================================
void clearMagnitudeColumns(ProcessConfig::RuntimeOptions& options) {
    options.extcat_mag_g_column_one_based = 0;
    options.extcat_mag_r_column_one_based = 0;
    options.extcat_mag_i_column_one_based = 0;
    options.extcat_mag_z_column_one_based = 0;
    options.extcat_mag_y_column_one_based = 0;
}

// ==========================================
// Function: Verify the legacy pass-through schema
// Method: Check the 18 + 1 + 29 contract and its terminal source field.
// ==========================================
void testLegacyLayout(int& failures) {
    ProcessConfig::RuntimeOptions options;
    PipelineCatalog::CatalogLayout layout;
    std::string error;
    expect(PipelineCatalog::resolveCatalogLayout(options, layout, error),
           "legacy layout should resolve: " + error, failures);
    expect(layout.external_columns == 18, "legacy external width should be 18",
           failures);
    expect(layout.ccd == 18, "legacy CCD index should be 18", failures);
    expect(layout.source_base == 19, "legacy source base should be 19",
           failures);
    expect(layout.source_columns == 29,
           "process_main suffix should contain 29 fields", failures);
    expect(layout.all_columns == 48, "legacy row width should be 48", failures);
    expect(layout.source.chi2 == 47, "legacy chi2 index should be 47", failures);
    expect(PipelineCatalog::describeCatalogLayout(layout).find("all_columns=48")
               != std::string::npos,
           "layout description should contain the total width", failures);
}

// ==========================================
// Function: Verify a compact pass-through schema with every magnitude present
// Method: Place all eight named fields in an eight-column raw catalog.
// ==========================================
void testCompactPassThroughLayout(int& failures) {
    ProcessConfig::RuntimeOptions options;
    options.extcat_total_columns = 8;
    options.extcat_ra_column_one_based = 1;
    options.extcat_dec_column_one_based = 2;
    options.extcat_mag_g_column_one_based = 3;
    options.extcat_mag_r_column_one_based = 4;
    options.extcat_mag_i_column_one_based = 5;
    options.extcat_mag_z_column_one_based = 6;
    options.extcat_mag_y_column_one_based = 7;
    options.extcat_zp_column_one_based = 8;

    PipelineCatalog::CatalogLayout layout;
    std::string error;
    expect(PipelineCatalog::resolveCatalogLayout(options, layout, error),
           "compact pass-through layout should resolve: " + error, failures);
    expect(layout.external_columns == 8,
           "compact pass-through external width should be 8", failures);
    expect(layout.ccd == 8 && layout.source_base == 9,
           "compact pass-through CCD/source offsets should be 8/9", failures);
    expect(layout.all_columns == 38,
           "compact pass-through row width should be 38", failures);
    expect(layout.external.ra == 0 && layout.external.dec == 1
               && layout.external.mag_g == std::size_t{2}
               && layout.external.mag_r == std::size_t{3}
               && layout.external.mag_i == std::size_t{4}
               && layout.external.mag_z == std::size_t{5}
               && layout.external.mag_y == std::size_t{6}
               && layout.external.zp == 7,
           "compact pass-through fields should map to indices 0 through 7",
           failures);
}

// ==========================================
// Function: Verify the minimal mandatory pass-through schema
// Method: Keep only RA/Dec/ZP named in three physical columns and mark every
//         magnitude absent through the raw zero sentinel.
// ==========================================
void testRequiredOnlyPassThroughLayout(int& failures) {
    ProcessConfig::RuntimeOptions options;
    options.extcat_total_columns = 3;
    options.extcat_ra_column_one_based = 1;
    options.extcat_dec_column_one_based = 2;
    options.extcat_zp_column_one_based = 3;
    clearMagnitudeColumns(options);

    PipelineCatalog::CatalogLayout layout;
    std::string error;
    expect(PipelineCatalog::resolveCatalogLayout(options, layout, error),
           "RA/Dec/ZP-only pass-through layout should resolve: " + error,
           failures);
    expect(layout.external_columns == 3 && layout.ccd == 3
               && layout.source_base == 4 && layout.all_columns == 33,
           "required-only layout should derive 3 + 1 + 29 columns", failures);
    expect(!layout.external.mag_g.has_value()
               && !layout.external.mag_r.has_value()
               && !layout.external.mag_i.has_value()
               && !layout.external.mag_z.has_value()
               && !layout.external.mag_y.has_value(),
           "all five magnitudes should be absent", failures);
    expect(PipelineCatalog::describeCatalogLayout(layout).find(
               "external_mag_i=none") != std::string::npos,
           "layout description should print none for an absent magnitude",
           failures);
}

// ==========================================
// Function: Verify ordered required-field projection
// Method: Project only the eight named fields in canonical order and check offsets.
// ==========================================
void testOrderedProjection(int& failures) {
    ProcessConfig::RuntimeOptions options = projectionOptions(
        {5, 6, 7, 9, 11, 13, 15, 17});
    PipelineCatalog::CatalogLayout layout;
    std::string error;
    expect(PipelineCatalog::resolveCatalogLayout(options, layout, error),
           "ordered projection should resolve: " + error, failures);
    expect(layout.external_columns == 8,
           "ordered projection external width should be 8", failures);
    expect(layout.ccd == 8, "ordered projection CCD index should be 8",
           failures);
    expect(layout.source_base == 9,
           "ordered projection source base should be 9", failures);
    expect(layout.all_columns == 38,
           "ordered projection row width should be 38", failures);
    expect(layout.external.ra == 0 && layout.external.dec == 1
               && layout.external.zp == 7,
           "ordered projection should map RA/Dec/ZP to 0/1/7", failures);
}

// ==========================================
// Function: Verify reordered projection and reader consumption
// Method: Move ZP/RA/Dec to the first columns and parse their effective positions.
// ==========================================
void testReorderedProjection(int& failures) {
    ProcessConfig::RuntimeOptions options = projectionOptions(
        {17, 5, 6, 7, 9, 11, 13, 15});
    PipelineCatalog::CatalogLayout layout;
    std::string error;
    expect(PipelineCatalog::resolveCatalogLayout(options, layout, error),
           "reordered projection should resolve: " + error, failures);
    expect(layout.external.zp == 0 && layout.external.ra == 1
               && layout.external.dec == 2,
           "reordered projection should map ZP/RA/Dec to 0/1/2", failures);
    expect(layout.external.mag_g == std::size_t{3}
               && layout.external.mag_r == std::size_t{4}
               && layout.external.mag_i == std::size_t{5}
               && layout.external.mag_z == std::size_t{6}
               && layout.external.mag_y == std::size_t{7},
           "reordered projection should map FD magnitudes to 3/4/5/6/7",
           failures);
    expect(ExternalCatalogReader::configure(layout, error),
           "reader should consume reordered layout: " + error, failures);

    ExternalCatalogReader::Record record;
    expect(ExternalCatalogReader::parseRecord("0.75 123.25 -45.5", record),
           "reader should parse the reordered leading fields", failures);
    expect(std::fabs(record.zp - 0.75) < 1.0e-12
               && std::fabs(record.ra - 123.25) < 1.0e-12
               && std::fabs(record.dec + 45.5) < 1.0e-12,
           "reader should assign reordered ZP/RA/Dec values correctly", failures);
}

// ==========================================
// Function: Verify extra physical columns remain unmodeled
// Method: Retain unrelated projected fields while resolving only eight named positions.
// ==========================================
void testExtraUnmodeledColumns(int& failures) {
    ProcessConfig::RuntimeOptions options = projectionOptions(
        {17, 5, 6, 1, 7, 8, 9, 11, 13, 15});
    PipelineCatalog::CatalogLayout layout;
    std::string error;
    expect(PipelineCatalog::resolveCatalogLayout(options, layout, error),
           "projection with extra unmodeled columns should resolve: " + error,
           failures);
    expect(layout.external_columns == 10 && layout.all_columns == 40,
           "extra fields should affect width but require no layout members",
           failures);
    expect(layout.external.zp == 0 && layout.external.ra == 1
               && layout.external.dec == 2
               && layout.external.mag_g == std::size_t{4}
               && layout.external.mag_r == std::size_t{6}
               && layout.external.mag_i == std::size_t{7}
               && layout.external.mag_z == std::size_t{8}
               && layout.external.mag_y == std::size_t{9},
           "named fields should skip over unmodeled projected columns", failures);
}

// ==========================================
// Function: Verify required and optional projection rules
// Method: Permit missing magnitudes while rejecting missing mandatory fields,
//         duplicate raw identities, and overlap among present named fields.
// ==========================================
void testProjectionRules(int& failures) {
    PipelineCatalog::CatalogLayout layout;
    std::string error;

    ProcessConfig::RuntimeOptions missing_ra = projectionOptions(
        {6, 7, 9, 11, 13, 15, 17});
    expect(!PipelineCatalog::resolveCatalogLayout(missing_ra, layout, error),
           "projection missing RA should fail", failures);

    ProcessConfig::RuntimeOptions missing_zp = projectionOptions(
        {5, 6, 7, 9, 11, 13, 15});
    expect(!PipelineCatalog::resolveCatalogLayout(missing_zp, layout, error),
           "projection missing ZP should fail", failures);

    ProcessConfig::RuntimeOptions missing_mag_i = projectionOptions(
        {17, 5, 6, 13});
    expect(PipelineCatalog::resolveCatalogLayout(missing_mag_i, layout, error),
           "projection missing optional mag_i should resolve: " + error,
           failures);
    expect(!layout.external.mag_i.has_value()
               && layout.external.mag_z == std::size_t{3},
           "projection omission should mark mag_i absent and retain mag_z",
           failures);

    ProcessConfig::RuntimeOptions duplicate = projectionOptions(
        {5, 6, 7, 9, 11, 11, 13, 15, 17});
    expect(!PipelineCatalog::resolveCatalogLayout(duplicate, layout, error),
           "projection with duplicate raw columns should fail", failures);

    ProcessConfig::RuntimeOptions overlap;
    overlap.extcat_mag_g_column_one_based =
        overlap.extcat_ra_column_one_based;
    expect(!PipelineCatalog::resolveCatalogLayout(overlap, layout, error),
           "two required meanings assigned to one raw column should fail",
           failures);
}

// ==========================================
// Function: Verify deterministic FD magnitude fallback
// Method: Remove available bands one at a time and require i-z-r-g-y priority,
//         followed by an explicit failure when no magnitude remains.
// ==========================================
void testFdMagnitudeSelection(int& failures) {
    ProcessConfig::RuntimeOptions options;
    PipelineCatalog::CatalogLayout layout;
    ProcessFD::MagnitudeSelection selection;
    std::string error;

    const auto expect_band = [&](const std::string& expected) {
        expect(PipelineCatalog::resolveCatalogLayout(options, layout, error),
               "magnitude-selection layout should resolve: " + error,
               failures);
        expect(ProcessFD::selectMagnitude(layout, selection, error),
               "FD magnitude selection should succeed: " + error, failures);
        expect(selection.band == expected,
               "FD should select " + expected + " band", failures);
    };

    expect_band("i");
    options.extcat_mag_i_column_one_based = 0;
    expect_band("z");
    options.extcat_mag_z_column_one_based = 0;
    expect_band("r");
    options.extcat_mag_r_column_one_based = 0;
    expect_band("g");
    options.extcat_mag_g_column_one_based = 0;
    expect_band("y");
    options.extcat_mag_y_column_one_based = 0;
    expect(PipelineCatalog::resolveCatalogLayout(options, layout, error),
           "layout without magnitudes should remain valid: " + error,
           failures);
    expect(!ProcessFD::selectMagnitude(layout, selection, error),
           "FD magnitude selection should fail when every band is absent",
           failures);
    expect(error.find("requires at least one external magnitude")
               != std::string::npos,
           "no-magnitude failure should explain the FD requirement", failures);
}

// ==========================================
// Function: Verify the selected fallback band drives FD catalog ingestion
// Method: Omit i, retain an out-of-range g and valid z, select z by priority,
//         then require the reader to accept and store the z magnitude.
// ==========================================
void testSelectedMagnitudeReader(int& failures) {
    ProcessConfig::RuntimeOptions options = projectionOptions(
        {17, 5, 6, 7, 13});
    PipelineCatalog::CatalogLayout layout;
    ProcessFD::MagnitudeSelection selection;
    std::string error;
    expect(PipelineCatalog::resolveCatalogLayout(options, layout, error),
           "selected-magnitude reader layout should resolve: " + error,
           failures);
    expect(ProcessFD::selectMagnitude(layout, selection, error)
               && selection.band == "z",
           "reader test should select z after missing i", failures);

    std::vector<float> row(layout.all_columns, 0.0f);
    row[layout.external.ra] = 10.0f;
    row[layout.external.dec] = 10.0f;
    row[layout.external.zp] = 0.5f;
    row[layout.external.mag_g.value()] = 5.0f;
    row[layout.external.mag_z.value()] = 22.5f;
    row[layout.ccd] = 1.0f;
    row[layout.source.pixx] = 1000.0f;
    row[layout.source.pixy] = 1000.0f;
    row[layout.source.h_flux] = 100.0f;
    row[layout.source.h_area] = 4.0f;
    row[layout.source.star] = 20.0f;
    row[layout.source.psf] = 1.0f;
    row[layout.source.snr_f] = 10.0f;
    row[layout.source.ra] = 10.0f;
    row[layout.source.dec] = 10.0f;
    row[layout.source.g1] = 0.1f;
    row[layout.source.g2] = 0.1f;
    row[layout.source.de] = 1.0f;
    row[layout.source.imax] = 10.0f;
    row[layout.source.jmax] = 10.0f;
    row[layout.source.flag] = 1.0f;

    const std::string filename =
        "tests/catalog_layout_selected_magnitude.tmp";
    {
        std::ofstream output(filename);
        output << "header\n" << numericRow(row) << '\n';
    }
    FDData data;
    data.reserve(1);
    ShearCatalogReader::readExposure(1, data, {filename}, layout,
                                     selection.column, 0);
    std::remove(filename.c_str());

    expect(data.ng == 1,
           "valid z magnitude should pass even when retained g is out of range",
           failures);
    expect(data.ng == 1 && std::fabs(data.star_mag[0] - 22.5f) < 1.0e-6f,
           "reader should store the selected z magnitude for star-bar use",
           failures);
}

// ==========================================
// Function: Verify exact runtime FD row width parsing
// Method: Accept N columns and reject both N-1 and N+1 rows.
// ==========================================
void testFdRowWidth(int& failures) {
    ProcessConfig::RuntimeOptions options = projectionOptions(
        {17, 5, 6, 7, 9, 11, 13, 15});
    PipelineCatalog::CatalogLayout layout;
    std::string error;
    expect(PipelineCatalog::resolveCatalogLayout(options, layout, error),
           "row-width test layout should resolve: " + error, failures);
    expect(layout.all_columns == 38,
           "eight-field FD row width should be 38", failures);

    std::vector<float> values;
    expect(ShearCatalogReader::parseCatalogRow(
               numericRow(layout.all_columns), layout.all_columns, values),
           "FD reader should accept the exact runtime row width", failures);
    expect(values.size() == layout.all_columns,
           "accepted FD row should retain every column", failures);
    expect(!ShearCatalogReader::parseCatalogRow(
               numericRow(layout.all_columns - 1), layout.all_columns, values),
           "FD reader should reject a short row", failures);
    expect(!ShearCatalogReader::parseCatalogRow(
               numericRow(layout.all_columns + 1), layout.all_columns, values),
           "FD reader should reject a row with trailing columns", failures);
}

}  // namespace

// ==========================================
// Function: Run CatalogLayout regression coverage
// Method: Execute legacy, projection, invalid-schema, reader, and row-width checks.
// ==========================================
int main() {
    int failures = 0;
    testLegacyLayout(failures);
    testCompactPassThroughLayout(failures);
    testRequiredOnlyPassThroughLayout(failures);
    testOrderedProjection(failures);
    testReorderedProjection(failures);
    testExtraUnmodeledColumns(failures);
    testProjectionRules(failures);
    testFdMagnitudeSelection(failures);
    testSelectedMagnitudeReader(failures);
    testFdRowWidth(failures);

    if (failures != 0) {
        std::cerr << failures << " catalog-layout test(s) failed\n";
        return 1;
    }
    std::cout << "CatalogLayout tests passed\n";
    return 0;
}
