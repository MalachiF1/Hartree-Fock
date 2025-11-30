# Makefile for the C++ Hartree-Fock Program

# Compiler
CXX = g++

MODE ?= release

PREPROCESSOR_DEFINES =

# Define the library name and potential installation paths for OpenBLAS
OPENBLAS_LIB_NAME = openblas
POSSIBLE_LIB_PATHS = /usr/lib /usr/local/lib /opt/OpenBLAS/lib

# Check for the OpenBLAS library file in the possible paths
# Note: The exact filename might vary (e.g., libopenblas.so.0, libopenblas.a)
FOUND_OPENBLAS_LIB := $(strip $(shell find $(POSSIBLE_LIB_PATHS) -name "lib$(OPENBLAS_LIB_NAME).so*" -print -quit))

# Conditional flags
ifeq ($(FOUND_OPENBLAS_LIB),)
    # OpenBLAS not found, no special flags needed, potentially link with default BLAS (e.g. standard libblas) or no BLAS
    BLAS_LIBS =
    BLAS_INCLUDE =

    $(warning OpenBLAS not found, building without specific OpenBLAS linkage)
else
    # OpenBLAS found, extract the path and set flags
    BLAS_LIB_PATH := $(dir $(FOUND_OPENBLAS_LIB))
    # Remove trailing slash from path for -L flag
    BLAS_LIB_DIR := $(patsubst %/,%,$(BLAS_LIB_PATH))
    # Set the linker and include flags
    BLAS_LIBS = -L$(BLAS_LIB_DIR) -l$(OPENBLAS_LIB_NAME)
    # Assuming standard include path relative to lib path, adjust if necessary
    BLAS_INCLUDE = -I$(BLAS_LIB_DIR)/../include
	# Also set preprocessor definitions to enable BLAS/LAPACK in Eigen
	PREPROCESSOR_DEFINES += -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACK
endif

ifeq ($(MODE), debug)
# Compiler flags
	CXXFLAGS = -std=c++20 $(PREPROCESSOR_DEFINES) -Iinclude -Ieigen -Inlohmann -Ifmt/include $(BLAS_INCLUDE) -fopenmp -pthread -fno-omit-frame-pointer -g -Wall -Wextra -Wpedantic -O0 -fno-inline -fsanitize=address
	LDFLAGS = -fopenmp $(BLAS_LIBS) -g
else ifeq ($(MODE), profile)
	CXXFLAGS = -std=c++20 $(PREPROCESSOR_DEFINES) -Iinclude -Ieigen -Inlohmann -Ifmt/include $(BLAS_INCLUDE) -fopenmp -pthread -fno-omit-frame-pointer -g -Wall -Wextra -Wpedantic -O3 -march=native -DNDEBUG
	LDFLAGS = -fopenmp $(BLAS_LIBS) -g
else ifeq ($(MODE), release)
# Linker flags
	CXXFLAGS = -std=c++20 $(PREPROCESSOR_DEFINES) -Iinclude -Ieigen -Inlohmann -Ifmt/include $(BLAS_INCLUDE) -fopenmp -g0 -DNDEBUG -march=native -O3 -w  -flto -funroll-loops -ffp-contract=fast -fno-math-errno -freciprocal-math -fno-trapping-math
 	LDFLAGS = -fopenmp $(BLAS_LIBS) -flto
else
	$(error "Invalid MODE specified. Use 'debug' or 'release'.")
endif

# Directories
SRCDIR = src
INCDIR = include
OBJDIR = obj
BINDIR = bin

# We need to build the fmt library as well (not header-only).
FMT_SRC = fmt/src/format.cc
FMT_OBJ = $(OBJDIR)/fmt.o

# Executable name
TARGET = $(BINDIR)/Hartree-Fock

# Automatically find all .cpp files in the source directory
SOURCES = $(wildcard $(SRCDIR)/*.cpp $(SRCDIR)/*/*.cpp)

# Generate corresponding object file names in the object directory
OBJECTS = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SOURCES))

# Default target: build the executable
all: $(TARGET)

# Rule to link the final executable
# The pipe '|' indicates an order-only prerequisite, ensuring the directory exists before linking.
$(TARGET): $(OBJECTS) $(FMT_OBJ) | $(BINDIR)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "Build complete. Run the program with: make run"

# Rule to compile the fmt library
$(FMT_OBJ): $(FMT_SRC)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Rule to compile a .cpp source file into a .o object file
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Create the output directories if they don't exist
$(BINDIR):
	mkdir -p $@

# Run the program after ensuring it's built
run: all
	./$(TARGET) $(filter-out run,$(MAKECMDGOALS))

profile: all
	perf record -a -g -e cycles -- ./$(TARGET) $(filter-out profile,$(MAKECMDGOALS))
	perf annotate -i perf.data --stdio "SCF::buildFockMatrix() [clone ._omp_fn.0]" > perf_annotate_buildFockMatrix.txt
	perf annotate -i perf.data --stdio "IntegralEngine::electronRepulsion(Basis const&, Shell const&, Shell const&, Shell const&, Shell const&, ElectronRepulsionTensor&)" > perf_annotate_electronRepulsion.txt
	perf annotate -i perf.data --stdio "Boys::calculateBoys(unsigned int, double, std::span<double, 18446744073709551615ul>)" > perf_annotate_boys.txt
	perf report --stdio > perf_report.txt
	@echo "Profiling complete. See perf_report.txt or use 'perf report' to view results."
	@echo "See perf_annotate_buildFockMatrix.txt, perf_annotate_electronRepulsion.txt, and perf_annotate_boys.txt for annotated functions."

# Clean up build artifacts
clean:
	rm -rf $(OBJDIR) $(BINDIR)
	@echo "Cleaned up build files."

# Phony targets are not real files
.PHONY: all clean run
