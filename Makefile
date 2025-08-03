# Makefile for the C++ Hartree-Fock Program

# Compiler
CXX = g++

# Debugging / optimization flags
# -Wall -Wextra / -O3 -march=native -flto (also add -flto in LDFLAGS)

# Compiler flags
CXXFLAGS = -std=c++17 -Iinclude -Ieigen -Inlohmann -fopenmp -O3 -march=native -pg

# Linker flags (e.g., for linking external libraries)
# link OpenBLAS statically
LDFLAGS = -Wl,-Bstatic -lopenblas -Wl,-Bdynamic

# Directories
SRCDIR = src
INCDIR = include
OBJDIR = obj
BINDIR = bin

# Executable name
TARGET = $(BINDIR)/hf_program

# Automatically find all .cpp files in the source directory
SOURCES = $(wildcard $(SRCDIR)/*.cpp)

# Generate corresponding object file names in the object directory
OBJECTS = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SOURCES))

# Default target: build the executable
all: $(TARGET)

# Rule to link the final executable
# The pipe '|' indicates an order-only prerequisite, ensuring the directory exists before linking.
$(TARGET): $(OBJECTS) | $(BINDIR)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "Build complete. Run the program with: make run"

# Rule to compile a .cpp source file into a .o object file
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Create the output directories if they don't exist
$(OBJDIR) $(BINDIR):
	mkdir -p $@

# Run the program after ensuring it's built
run: all
	./$(TARGET)

# Clean up build artifacts
clean:
	rm -rf $(OBJDIR) $(BINDIR)
	@echo "Cleaned up build files."

# Phony targets are not real files
.PHONY: all clean run

