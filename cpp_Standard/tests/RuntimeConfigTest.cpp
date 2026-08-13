#include "CatalogLayout.hpp"
#include "RuntimeConfig.hpp"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

// ==========================================
// Function: expect
// Method: Accumulate focused runtime-config test failures while emitting only
//         violated invariants for readable command-line regression output.
// ==========================================
void expect(bool condition, const std::string& message, int& failures) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

// ==========================================
// Function: parseArguments
// Method: Convert mutable string storage to argv form and exercise the exact
//         production phase-B CLI parser without a subprocess.
// ==========================================
bool parseArguments(const std::vector<std::string>& arguments,
                    RuntimeConfig& config, std::string& error) {
    std::vector<std::string> storage = arguments;
    std::vector<char*> argv;
    argv.reserve(storage.size());
    for (std::string& argument : storage) {
        argv.push_back(argument.data());
    }
    return parseRuntimeCommandLine(static_cast<int>(argv.size()), argv.data(),
                                   config, error);
}

// ==========================================
// Function: testDefaults
// Method: Verify the effective fallback value set, including the user's
//         disabled legacy path limit and all runtime geometry defaults.
// ==========================================
void testDefaults(int& failures) {
    const RuntimeConfig config = makeDefaultRuntimeConfig();
    expect(!config.process.run_process_extcat
               && config.process.run_process_init
               && config.process.run_process_main
               && config.process.run_process_rearr
               && config.process.run_process_fd,
           "compiled Standard phase defaults should be preserved", failures);
    expect(config.init.f77_max_path == 0,
           "user-configured F77 path default should remain disabled", failures);
    expect(config.lensing.process_stage == 223092870,
           "all nine numerical stages should remain enabled by default", failures);
    expect(config.lensing.nmax_chip == 62 && config.lensing.chipnx == 2046
               && config.lensing.chipny == 4094,
           "camera defaults should seed runtime geometry", failures);
}

// ==========================================
// Function: testIniSectionsAndTypes
// Method: Exercise every INI value family, quoted strings, aliases, comments,
//         lists, datasets, optional magnitude sentinels, and schema derivation.
// ==========================================
void testIniSectionsAndTypes(int& failures) {
    RuntimeConfig config = makeDefaultRuntimeConfig();
    const std::string text = R"ini(
# full parser coverage
[process]
run_process_extcat = yes
run_process_init = false
run_process_main = on
run_process_rearr = 1
run_process_fd = 0
expo_list = "/data/catalog path/expo.list" ; quoted path
rearr_output_directory = selected
rearr_output_base_directory = /tmp/rearr
rearranged_expo_list_filename = custom.list
rearranged_expo_list_directory = /tmp/lists
fd_expo_list = /tmp/fd.list
fd_output_directory = fd-selected
fd_output_base_directory = /tmp/fd

[extcat]
input_directory = /data/raw
output_directory = /data/tiles
filename_tokens = "tile,part", source_
recursive = no
delimiter = comma
header_mode = present
malformed_policy = skip
existing_policy = overwrite
chunk_mib = 128
total_columns = 20
use_explicit_columns = true
input_columns = 5,6,17,7,9
use_explicit_coordinate_columns = true
ra_column = 5
dec_column = 6
mag_g_column = 7
mag_r_column = 9
mag_i_column = 0
mag_z_column = 0
mag_y_column = 0
zp_column = 17

[init]
science_root = /archive/science
dq_root = /archive/dq
output_root = /work/output
datasets = gband:c4d_, rband:r4d_
contains = v1, "release,2"
existing = resume
f77_max_path = 210

[lensing]
astrometry_trivial = 1
process_stage = 19
include_flat = 1
include_mask = 3
astrometry_cat = /catalog/gaia
source_cat = "/catalog/source tiles"
flat_path = /calibration/flat
psf_path = /calibration/psf
ext_cat = 1
ext_psf = 1
ccd_split = 1
deblending = 0
psf_type = 2
psf_ms = 1
gal_smooth = 3
star_smooth = 4
pixel_size = 0.31
nmax_chip = 71
chipnx = 2500
chipny = 4500
)ini";
    std::string error;
    expect(applyRuntimeConfigText(text, "test.ini", config, error),
           "valid full INI should parse: " + error, failures);
    expect(config.process.run_process_extcat
               && !config.process.run_process_init
               && config.process.run_process_main
               && config.process.run_process_rearr
               && !config.process.run_process_fd,
           "INI booleans should accept yes/no/on/numeric forms", failures);
    expect(config.process.expo_list == "/data/catalog path/expo.list",
           "quoted path and inline comment should decode", failures);
    expect(config.extcat.filename_tokens.size() == 2
               && config.extcat.filename_tokens[0] == "tile,part"
               && config.extcat.filename_tokens[1] == "source_",
           "quoted comma list values should remain one item", failures);
    expect(config.extcat.input_columns_one_based
               == std::vector<std::size_t>({5, 6, 17, 7, 9}),
           "explicit projection should preserve configured order", failures);
    expect(config.init.datasets.size() == 2
               && config.init.datasets[1].target == "rband"
               && config.init.contains[1] == "release,2",
           "dataset and token lists should parse", failures);
    expect(config.extcat.output_directory == "/catalog/source tiles",
           "lensing.source_cat should alias the authoritative extcat output", failures);
    expect(config.lensing.nmax_chip == 71 && config.lensing.chipnx == 2500
               && config.lensing.chipny == 4500
               && std::fabs(config.lensing.pixel_size - 0.31) < 1.0e-12,
           "runtime geometry and pixel scale should parse", failures);

    PipelineCatalog::CatalogLayout layout;
    expect(PipelineCatalog::resolveCatalogLayout(config, layout, error),
           "projection layout should resolve after INI parsing: " + error,
           failures);
    expect(layout.external_columns == 5 && layout.external.ra == 0
               && layout.external.dec == 1 && layout.external.zp == 2
               && layout.external.mag_g == std::size_t{3}
               && layout.external.mag_r == std::size_t{4},
           "layout should consume the nested extcat runtime section", failures);
}

// ==========================================
// Function: testTransactionalErrors
// Method: Ensure unknown keys, malformed types, unterminated quotes, and keys
//         before sections fail with source context and do not partially commit.
// ==========================================
void testTransactionalErrors(int& failures) {
    const std::vector<std::string> invalid_texts = {
        "[lensing]\nchipnx=1234\nunknown=1\n",
        "[lensing]\nchipnx=abc\n",
        "[process]\nexpo_list=\"unterminated\n",
        "chipnx=10\n",
        "[unknown]\nvalue=1\n",
    };
    for (const std::string& text : invalid_texts) {
        RuntimeConfig config = makeDefaultRuntimeConfig();
        std::string error;
        expect(!applyRuntimeConfigText(text, "invalid.ini", config, error),
               "invalid INI should fail", failures);
        expect(error.find("file: invalid.ini") != std::string::npos
                   && error.find("reason:") != std::string::npos,
               "invalid INI should carry structured source diagnostics", failures);
        expect(config.lensing.chipnx == 2046,
               "failed INI should not commit earlier values", failures);
    }
}

// ==========================================
// Function: testPrecedenceAndConfigScan
// Method: Demonstrate default < INI < CLI, both CLI syntaxes, list replacement,
//         positional compatibility, and phase-A last-config-path behavior.
// ==========================================
void testPrecedenceAndConfigScan(int& failures) {
    RuntimeConfig config = makeDefaultRuntimeConfig();
    std::string error;
    expect(applyRuntimeConfigText(
               "[process]\nrun_process_main=false\nexpo_list=file.list\n"
               "[init]\ndatasets=file:pfx\ncontains=filetoken\n",
               "precedence.ini", config, error),
           "precedence INI should parse: " + error, failures);
    expect(parseArguments(
               {"fq", "--config", "precedence.ini", "--run-main=true",
                "--dataset", "cli:prefix", "--contains=cli-token",
                "--expo-list", "cli.list"},
               config, error),
           "CLI overrides should parse: " + error, failures);
    expect(config.process.run_process_main
               && config.process.expo_list == "cli.list",
           "CLI scalar values should override INI values", failures);
    expect(config.init.datasets.size() == 1
               && config.init.datasets[0].target == "cli"
               && config.init.contains == std::vector<std::string>({"cli-token"}),
           "first CLI list option should replace the INI list", failures);

    std::vector<std::string> storage = {
        "fq", "--config=first.ini", "--help", "--config", "second.ini"};
    std::vector<char*> argv;
    for (std::string& argument : storage) argv.push_back(argument.data());
    std::string path;
    bool help = false;
    expect(findRuntimeConfigPath(static_cast<int>(argv.size()), argv.data(),
                                 path, help, error),
           "phase-A scan should succeed: " + error, failures);
    expect(path == "second.ini" && help,
           "last --config should win and --help should be detected", failures);
}

// ==========================================
// Function: testValidationAndSchemas
// Method: Cover positive-only geometry without an instrument ceiling, enum and
//         dependency errors, FD external-catalog gating, and internal rearr layout.
// ==========================================
void testValidationAndSchemas(int& failures) {
    RuntimeConfig config = makeDefaultRuntimeConfig();
    std::string error;
    expect(validateRuntimeConfig(config, error),
           "compiled defaults should validate: " + error, failures);

    config.process.run_process_extcat = false;
    config.process.run_process_init = false;
    config.process.run_process_main = false;
    config.process.run_process_rearr = false;
    config.process.run_process_fd = false;
    expect(validateRuntimeConfig(config, error),
           "an all-disabled configuration should support parse-only dry runs: "
               + error,
           failures);
    config = makeDefaultRuntimeConfig();

    config.lensing.nmax_chip = 1000;
    config.lensing.chipnx = 7;
    config.lensing.chipny = 9;
    expect(validateRuntimeConfig(config, error),
           "positive geometry should have no instrument hard limit: " + error,
           failures);

    config.lensing.ext_cat = 0;
    config.process.run_process_fd = false;
    config.extcat.output_directory.clear();
    expect(validateRuntimeConfig(config, error),
           "internal main/rearr mode should not require external catalog paths: "
               + error,
           failures);
    PipelineCatalog::RearrCatalogSchema schema;
    expect(PipelineCatalog::resolveRearrCatalogSchema(config, nullptr, schema,
                                                       error),
           "internal rearr schema should resolve without external layout: " + error,
           failures);
    expect(schema.all_columns == 30 && schema.ra_column == 13
               && schema.dec_column == 14,
           "internal rearr schema should be CCD + 29 source fields", failures);

    config.process.run_process_fd = true;
    expect(!validateRuntimeConfig(config, error)
               && error.find("requires ext_cat=1") != std::string::npos,
           "FD should fail early for ext_cat=0", failures);
    config.process.run_process_fd = false;
    config.lensing.process_stage = 23;
    expect(!validateRuntimeConfig(config, error)
               && error.find("requires Stage 8") != std::string::npos,
           "Stage 9 without Stage 8 should fail validation", failures);
}

// ==========================================
// Function: testWriteOnceStore
// Method: Verify deep consumers see the finalized value and a second
//         initialization cannot mutate it.
// ==========================================
void testWriteOnceStore(int& failures) {
    RuntimeConfig config = makeDefaultRuntimeConfig();
    config.lensing.chipnx = 3333;
    std::string error;
    expect(RuntimeConfigStore::initialize(config, error),
           "first RuntimeConfigStore initialization should succeed", failures);
    expect(RuntimeConfigStore::isInitialized()
               && RuntimeConfigStore::get().lensing.chipnx == 3333,
           "store should expose the finalized immutable value", failures);
    config.lensing.chipnx = 4444;
    expect(!RuntimeConfigStore::initialize(config, error)
               && RuntimeConfigStore::get().lensing.chipnx == 3333,
           "second store initialization should fail without mutation", failures);
}

}  // namespace

// ==========================================
// Function: main
// Method: Run parser, precedence, validation, schema, and store regressions in
//         deterministic order and return nonzero on any failed invariant.
// ==========================================
int main() {
    int failures = 0;
    testDefaults(failures);
    testIniSectionsAndTypes(failures);
    testTransactionalErrors(failures);
    testPrecedenceAndConfigScan(failures);
    testValidationAndSchemas(failures);
    testWriteOnceStore(failures);
    if (failures != 0) {
        std::cerr << failures << " RuntimeConfig test(s) failed\n";
        return 1;
    }
    std::cout << "RuntimeConfig tests passed\n";
    return 0;
}
