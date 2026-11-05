#include "CatalogLayout.hpp"
#include "RuntimeConfig.hpp"
#include "process_astrocat/process_astrocat.hpp"
#include "process_extcat/process_extcat.hpp"
#include "process_fd/process_fd.hpp"
#include "process_init/process_init.hpp"
#include "general/ExposureList.hpp"
#include "general/MPIScheduler.hpp"
#include "general/MPIUtils.hpp"
#include "general/NumericalRecipes.hpp"
#include "general/PathUtils.hpp"
#include "process_main/process_main.hpp"
#include "process_rearr/process_rearr.hpp"
#include <mpi.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace {
// All ranks follow the same startup/failure path, including INI I/O failures.
int allSucceeded(int local) {
    int global=0;
    MPI_Allreduce(&local,&global,1,MPI_INT,MPI_MIN,MPIScheduler::state.communicator);
    return global;
}
bool loadAndBroadcastConfigText(const std::string& path,int rank,
                                std::string& text,std::string& error) {
    int read_ok=1;
    if (rank==0) {
        std::ifstream input(path,std::ios::binary);
        if (!input.is_open()) { error="cannot open runtime config file: "+path; read_ok=0; }
        else {
            std::ostringstream buffer; buffer<<input.rdbuf();
            if (!input.good() && !input.eof()) { error="I/O failure while reading runtime config file: "+path; read_ok=0; }
            else text=buffer.str();
        }
    }
    MPI_Bcast(&read_ok,1,MPI_INT,0,MPIScheduler::state.communicator);
    if (!read_ok) return false;
    return MPIUtils::broadcastString(text,0,error);
}
void printUsage(const char* name) {
    std::cerr<<"Usage: "<<name<<" [--config PATH]\n"
             <<"No arguments: use compiled config/*.hpp defaults.\n"
             <<"--config PATH (or --config=PATH): override defaults with an INI file.\n"
             <<"No implicit pipeline.ini lookup. All other parameters belong in INI.\n";
}
std::string resolveExposureList(const RuntimeConfig& config,const InitConfig::DatasetSpec& dataset) {
    std::filesystem::path path;
    if (!config.process.expo_list.empty()) path=config.process.expo_list;
    else {
        if (config.init.output_root.empty())
            throw std::runtime_error("expo_list is absent and output_root cannot derive its default");
        path=std::filesystem::path(config.init.output_root)/("expo_"+dataset.target+".list");
    }
    return PathUtils::normalizedAbsolute(path).string();
}
std::string deriveDatasetRootFromExpoList(const std::string& list) {
    std::vector<ExposureList::Entry> exposures; std::string error;
    if (!ExposureList::loadPipelineList(list,exposures,0,error)) return "";
    for (const auto& exposure:exposures) {
        std::vector<std::string> images;
        if (!ExposureList::loadPathList(exposure.path,images,0,error)) continue;
        std::filesystem::path root;
        if (!PathUtils::parentAtLevel(images.front(),3,root,error)) continue;
        return PathUtils::normalizedAbsolute(root).string();
    }
    return "";
}
} // namespace

int main(int argc,char* argv[]) {
    MPIScheduler::init(argc,argv);
    const int rank=MPIScheduler::state.rank;
    if (rank==0) std::cout<<"MPI Init Done..."<<std::endl;
    std::string path,text,error; bool unused_help=false;
    int ok=allSucceeded(findRuntimeConfigPath(argc,argv,path,unused_help,error)?1:0);
    RuntimeConfig config=makeDefaultRuntimeConfig();
    if (ok && !path.empty()) {
        ok=allSucceeded(loadAndBroadcastConfigText(path,rank,text,error)?1:0);
        if (ok) ok=allSucceeded(applyRuntimeConfigText(text,path,config,error)?1:0);
    }
    if (ok) ok=allSucceeded(validateRuntimeConfig(config,error)?1:0);
    if (ok) ok=allSucceeded(RuntimeConfigStore::initialize(config,error)?1:0);
    if (!ok) {
        if (rank==0) {
            std::cerr<<"Config error:\n  "<<(error.empty()?"startup failed on another MPI rank":error)<<std::endl;
            printUsage(argv[0]);
        }
        MPIScheduler::finalize(); return 2;
    }
    const RuntimeConfig& runtime_config=RuntimeConfigStore::get();
    int return_code=0;

    std::optional<PipelineCatalog::CatalogLayout> external_layout;
    std::string schema_error;
    int schema_ok=1;
    if (runtime_config.lensing.ext_cat==1) {
        PipelineCatalog::CatalogLayout layout;
        schema_ok=allSucceeded(PipelineCatalog::resolveCatalogLayout(runtime_config,layout,schema_error)?1:0);
        if (schema_ok) external_layout=std::move(layout);
    }
    PipelineCatalog::RearrCatalogSchema rearr_schema;
    if (schema_ok && runtime_config.process.run_process_rearr) {
        schema_ok=allSucceeded(PipelineCatalog::resolveRearrCatalogSchema(runtime_config,
            external_layout?&*external_layout:nullptr,rearr_schema,schema_error)?1:0);
    }

    if (!schema_ok) {
        if (rank==0) std::cerr<<"Catalog schema error: "
            <<(schema_error.empty()?"resolution failed on another rank":schema_error)<<std::endl;
        MPIScheduler::finalize(); return 2;
    }
    if (rank==0 && external_layout) std::cout<<PipelineCatalog::describeCatalogLayout(*external_layout)<<std::endl;
    if (runtime_config.process.run_process_astrocat) {
        if (rank==0) std::cout<<"Running process_astrocat before process_extcat"<<std::endl;
        return_code=process_astrocat(runtime_config);
        if (!return_code) MPIScheduler::barrier();
    }
    if (!return_code && runtime_config.process.run_process_extcat) {
        if (rank==0) std::cout<<"Running process_extcat before all dataset phases"<<std::endl;
        return_code=process_extcat(runtime_config);
        if (!return_code) MPIScheduler::barrier();
    }
    bool rng_initialized=false;
    for (std::size_t index=0; index<runtime_config.init.datasets.size() && !return_code
         && (runtime_config.process.run_process_init || runtime_config.process.run_process_main
             || runtime_config.process.run_process_rearr || runtime_config.process.run_process_fd); ++index) {
        const auto& dataset=runtime_config.init.datasets[index];
        if (rank==0) std::cout<<"Dataset "<<(index+1)<<"/"<<runtime_config.init.datasets.size()
            <<": target="<<dataset.target<<" prefix="<<dataset.prefix<<std::endl;
        std::string generated_exposure_list;
        if (runtime_config.process.run_process_init)
            return_code=process_init(runtime_config,dataset,generated_exposure_list);
        std::string selected_exposure_list;
        if (!return_code && (runtime_config.process.run_process_main || runtime_config.process.run_process_rearr
                             || runtime_config.process.run_process_fd)) {
            if (runtime_config.process.run_process_init) {
                selected_exposure_list=generated_exposure_list;
                if (rank==0) std::cout<<"Downstream phases will use process_init output: "<<selected_exposure_list<<std::endl;
                MPIScheduler::barrier();
            } else {
                int local_ok=1; std::string path_error;
                try { selected_exposure_list=resolveExposureList(runtime_config,dataset); }
                catch (const std::exception& ex) { local_ok=0; path_error=ex.what(); }
                if (!allSucceeded(local_ok)) {
                    if (rank==0) std::cerr<<"Exposure-list error: "
                        <<(path_error.empty()?"resolution failed on another rank":path_error)<<std::endl;
                    return_code=2;
                }
            }
        }
        if (!return_code && runtime_config.process.run_process_main) {
            if (!rng_initialized) {
                const unsigned int seed=NumericalRecipes::initializeRan1Seed(rank,MPIScheduler::state.size);
                std::cout<<"RNG_SEED rank seed: "<<rank<<" "<<seed<<std::endl;
                MPIScheduler::barrier(); rng_initialized=true;
            }
            return_code=process_main(selected_exposure_list,runtime_config,external_layout?&*external_layout:nullptr);
        }
        if (!return_code && runtime_config.process.run_process_rearr) {
            MPIScheduler::barrier();
            if (rank==0) std::cout<<"Running process_rearr"<<std::endl;
            return_code=process_rearr(selected_exposure_list,runtime_config,rearr_schema);
        }
        if (!return_code && runtime_config.process.run_process_fd) {
            MPIScheduler::barrier();
            if (!rng_initialized) {
                NumericalRecipes::initializeRan1Seed(rank,MPIScheduler::state.size);
                MPIScheduler::barrier(); rng_initialized=true;
            }
            std::string fd_expo_list;
            if (runtime_config.process.fd_expo_list.empty()) {
                const std::filesystem::path list_dir=runtime_config.process.rearranged_expo_list_directory.empty()
                    ?std::filesystem::path(selected_exposure_list).parent_path()
                    :std::filesystem::path(runtime_config.process.rearranged_expo_list_directory);
                fd_expo_list=std::filesystem::absolute(list_dir/runtime_config.process.rearranged_expo_list_filename)
                    .lexically_normal().string();
            } else fd_expo_list=runtime_config.process.fd_expo_list;
            const auto dataset_root=deriveDatasetRootFromExpoList(selected_exposure_list);
            if (rank==0) std::cout<<"Running process_fd: expo_list="<<fd_expo_list<<" dataset_root="<<dataset_root<<std::endl;
            return_code=process_fd(fd_expo_list,runtime_config,dataset_root,*external_layout);
        }
        if (!return_code && index+1<runtime_config.init.datasets.size()) MPIScheduler::barrier();
    }
    MPIScheduler::finalize(); return return_code;
}
