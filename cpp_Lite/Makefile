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
       src/process_extcat/process_extcat.cpp \
       src/process_rearr/process_rearr.cpp \
       src/process_rearr/CatalogRearranger.cpp \
       src/process_fd/process_fd.cpp \
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
       src/process_main/ShearMeasurement.cpp \
       src/process_main/ExposureInfo.cpp \
       src/process_main/CatalogCombiner.cpp \
       src/process_main/FitsIO.cpp \
       src/process_main/LinearSolve.cpp \
       src/process_main/UniversalUtils.cpp \
       src/process_main/ImageProcessing.cpp \
       src/process_main/NumericalRecipes.cpp \
       src/process_main/MPIScheduler.cpp \
       src/process_main/ExStar.cpp

OBJS = $(SRCS:.cpp=.o)
DEPS = $(OBJS:.o=.d)
TARGET = Fourier_Quad_Pipe

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(DEPS) $(TARGET)

-include $(DEPS)

.PHONY: all clean
