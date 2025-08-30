# Makefile for the C++ Hartree-Fock Program

# Compiler
CXX = g++

# Debugging / optimization flags
# -g -Wall -Wextra / -O3 -march=native -flto (also add -flto in LDFLAGS)

# Compiler flags
CXXFLAGS = -std=c++20 -Iinclude -Ieigen -Inlohmann -Ifmt/include -fopenmp -pg -O3 -march=native -Wall -Wextra

# Linker flags (e.g., for linking external libraries)
# link OpenBLAS statically
LDFLAGS = -fopenmp -Wl,-Bstatic -lopenblas -Wl,-Bdynamic -pg

# Directories
SRCDIR = src
INCDIR = include
OBJDIR = obj
BINDIR = bin

# We need to build the fmt library as well (not header-only).
FMT_SRC = fmt/src/format.cc
FMT_OBJ = $(OBJDIR)/fmt.o

# Executable name
TARGET = $(BINDIR)/hf_program

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
	./$(TARGET)
	gprof -b bin/hf_program gmon.out > profiling_output.txt
	conda run -n sci gprof2dot profiling_output.txt | dot -Tpng -o callgraph.png
	@echo "Profiling complete. Check profiling_output.txt and callgraph.png for results."

# Clean up build artifacts
clean:
	rm -rf $(OBJDIR) $(BINDIR)
	@echo "Cleaned up build files."

# Phony targets are not real files
.PHONY: all clean run
