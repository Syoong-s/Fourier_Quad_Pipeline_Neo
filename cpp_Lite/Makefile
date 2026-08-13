CXX = mpicxx
STACK_PREFIX ?=
EIGEN_INCLUDE ?=

CPPFLAGS = -Iinclude -Iconfig -Iinclude/process_extcat -Iinclude/process_init -Iinclude/process_main -Iinclude/process_rearr -Iinclude/process_fd
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra -MMD -MP
LDFLAGS =
LDLIBS = -lcfitsio -lfftw3 -lfftw3f -llapack -lblas -lm

ifneq ($(strip $(STACK_PREFIX)),)
CPPFLAGS += -I$(STACK_PREFIX)/include -I$(STACK_PREFIX)/include/eigen3
LDFLAGS += -L$(STACK_PREFIX)/lib -Wl,-rpath,$(STACK_PREFIX)/lib
endif

ifneq ($(strip $(EIGEN_INCLUDE)),)
CPPFLAGS += -I$(EIGEN_INCLUDE)
endif

SRCS = main.cpp \
       src/RuntimeConfig.cpp \
       src/CatalogLayout.cpp \
       src/process_extcat/process_extcat.cpp \
       src/process_rearr/process_rearr.cpp \
       src/process_rearr/CatalogRearranger.cpp \
       src/process_fd/process_fd.cpp \
       src/process_fd/FDMagnitudeSelector.cpp \
       src/process_fd/ShearCatalogReader.cpp \
       src/process_fd/StarCutCalculator.cpp \
       src/process_fd/KMeansClusterer.cpp \
       src/process_fd/FDMeasurement.cpp \
       src/process_fd/QuadraticFitting.cpp \
       src/process_main/process_main.cpp \
       src/process_main/ExternalCatalogReader.cpp \
       src/process_init/process_init.cpp \
       src/process_init/Initializer.cpp \
       src/process_init/FitsExtractor.cpp \
       src/process_main/PreProcess.cpp \
       src/process_main/Astrometry.cpp \
       src/process_main/SourceExtractor.cpp \
       src/process_main/FourierTransformSt1.cpp \
       src/process_main/FourierTransformSt2.cpp \
       src/process_main/PSFModel.cpp \
       src/process_main/CurvatureSizeMeasurement.cpp \
       src/process_main/PointSourceStatistics.cpp \
       src/process_main/ShearMeasurement.cpp \
       src/process_main/ExposureInfo.cpp \
       src/process_main/CatalogCombiner.cpp \
       src/process_main/FitsIO.cpp \
       src/process_main/LinearSolve.cpp \
       src/process_main/Universalblock.cpp \
       src/process_main/UniversalUtils.cpp \
       src/process_main/ImageProcessing.cpp \
       src/process_main/NumericalRecipes.cpp \
       src/process_main/MPIScheduler.cpp \
       src/process_main/ExStar.cpp

OBJS = $(SRCS:.cpp=.o)
DEPS = $(OBJS:.o=.d)
TARGET = Fourier_Quad_Pipe
POINT_SOURCE_TEST = tests/test_point_source_statistics
STAMP_CUBE_IO_TEST = tests/test_stamp_cube_io
UNIVERSALBLOCK_TEST = tests/test_universalblock
CATALOG_LAYOUT_TEST = tests/test_catalog_layout
CATALOG_ROW_COUNT_TEST = tests/test_catalog_row_count
PSF_MODEL_STATE_TEST = tests/test_psf_model_state
RUNTIME_CONFIG_TEST = tests/test_runtime_config

all: $(TARGET)

test-point-source-statistics: $(POINT_SOURCE_TEST)
	./$(POINT_SOURCE_TEST)

test-stamp-cube-io: $(STAMP_CUBE_IO_TEST)
	./$(STAMP_CUBE_IO_TEST)

test-universalblock: $(UNIVERSALBLOCK_TEST)
	./$(UNIVERSALBLOCK_TEST)

test-catalog-layout: $(CATALOG_LAYOUT_TEST)
	./$(CATALOG_LAYOUT_TEST)

test-catalog-row-count: $(CATALOG_ROW_COUNT_TEST)
	./$(CATALOG_ROW_COUNT_TEST)

test-psf-model-state: $(PSF_MODEL_STATE_TEST)
	./$(PSF_MODEL_STATE_TEST)

test-runtime-config: $(RUNTIME_CONFIG_TEST)
	./$(RUNTIME_CONFIG_TEST)

$(POINT_SOURCE_TEST): tests/PointSourceStatisticsTest.cpp \
                     src/process_main/PointSourceStatistics.cpp \
                     include/process_main/PointSourceStatistics.hpp \
                     config/LensingConfig.hpp
	$(CXX) $(CPPFLAGS) $(filter-out -MMD -MP,$(CXXFLAGS)) \
		tests/PointSourceStatisticsTest.cpp src/process_main/PointSourceStatistics.cpp -o $@

$(STAMP_CUBE_IO_TEST): tests/StampCubeIOTest.cpp \
                      src/process_main/FitsIO.cpp \
                      include/process_main/FitsIO.hpp \
                      include/process_main/OutputFile.hpp
	$(CXX) $(CPPFLAGS) $(filter-out -MMD -MP,$(CXXFLAGS)) \
		tests/StampCubeIOTest.cpp src/process_main/FitsIO.cpp -o $@ -lcfitsio

$(UNIVERSALBLOCK_TEST): tests/UniversalblockTest.cpp \
                       src/process_main/Universalblock.cpp \
                       src/process_main/UniversalUtils.cpp \
                       src/process_main/FitsIO.cpp \
                       src/process_main/LinearSolve.cpp \
                       include/process_main/Universalblock.hpp \
                       include/process_main/FitsIO.hpp \
                       include/OutputLayout.hpp
	$(CXX) $(CPPFLAGS) $(filter-out -MMD -MP,$(CXXFLAGS)) \
		tests/UniversalblockTest.cpp \
		src/process_main/Universalblock.cpp \
		src/process_main/UniversalUtils.cpp \
		src/process_main/FitsIO.cpp \
		src/process_main/LinearSolve.cpp \
		$(LDFLAGS) -o $@ $(LDLIBS)

$(CATALOG_LAYOUT_TEST): tests/CatalogLayoutTest.cpp \
                       src/RuntimeConfig.cpp \
                       src/CatalogLayout.cpp \
                       src/process_main/ExternalCatalogReader.cpp \
                       src/process_fd/FDMagnitudeSelector.cpp \
                       src/process_fd/ShearCatalogReader.cpp \
                       include/CatalogLayout.hpp \
                       include/process_main/ExternalCatalogReader.hpp \
                       include/process_fd/FDMagnitudeSelector.hpp \
                       include/process_fd/ShearCatalogReader.hpp
	$(CXX) $(CPPFLAGS) $(filter-out -MMD -MP,$(CXXFLAGS)) \
		tests/CatalogLayoutTest.cpp \
		src/RuntimeConfig.cpp \
		src/CatalogLayout.cpp \
		src/process_main/ExternalCatalogReader.cpp \
		src/process_fd/FDMagnitudeSelector.cpp \
		src/process_fd/ShearCatalogReader.cpp -o $@

$(CATALOG_ROW_COUNT_TEST): tests/CatalogRowCountTest.cpp \
                          include/process_main/CatalogRowCount.hpp \
                          include/process_main/MPIFailure.hpp
	$(CXX) $(CPPFLAGS) $(filter-out -MMD -MP,$(CXXFLAGS)) \
		tests/CatalogRowCountTest.cpp -o $@

$(PSF_MODEL_STATE_TEST): tests/PSFModelStateTest.cpp \
                        include/process_main/PSFModelState.hpp \
                        config/LensingConfig.hpp
	$(CXX) $(CPPFLAGS) $(filter-out -MMD -MP,$(CXXFLAGS)) \
		tests/PSFModelStateTest.cpp -o $@

$(RUNTIME_CONFIG_TEST): tests/RuntimeConfigTest.cpp \
                       src/RuntimeConfig.cpp \
                       src/CatalogLayout.cpp \
                       include/RuntimeConfig.hpp \
                       include/CatalogLayout.hpp \
                       config/ProcessConfig.hpp \
                       config/ExtCatConfig.hpp \
                       config/InitConfig.hpp \
                       config/LensingConfig.hpp
	$(CXX) $(CPPFLAGS) $(filter-out -MMD -MP,$(CXXFLAGS)) \
		tests/RuntimeConfigTest.cpp src/RuntimeConfig.cpp \
		src/CatalogLayout.cpp -o $@

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(DEPS) $(TARGET) $(POINT_SOURCE_TEST) \
	      $(STAMP_CUBE_IO_TEST) $(UNIVERSALBLOCK_TEST) \
	      $(CATALOG_LAYOUT_TEST) $(CATALOG_ROW_COUNT_TEST) \
	      $(PSF_MODEL_STATE_TEST) $(RUNTIME_CONFIG_TEST) \
	      tests/catalog_layout_selected_magnitude.tmp

-include $(DEPS)

.PHONY: all clean test-point-source-statistics test-stamp-cube-io \
		test-universalblock test-catalog-layout test-catalog-row-count \
		test-psf-model-state test-runtime-config
