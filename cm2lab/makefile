# Modern Makefile for cm2lab
# Fortran computational physics code

# Compiler configuration
FC := gfortran
AR := ar
RANLIB := ranlib

# Directory structure
SRCDIR := .
BUILDDIR := build
TESTDIR := test

# Source files
SOURCES := cm2lab.f gauss.f90
OBJECTS := $(SOURCES:%.f=$(BUILDDIR)/%.o)
OBJECTS := $(OBJECTS:%.f90=$(BUILDDIR)/%.o)
TARGET := cm2lab

# Compiler flags
FFLAGS := -O2 -Wall -Wtabs -ffixed-form -fno-align-commons
FFLAGS_DEBUG := -g -fbounds-check -fbacktrace -Wall -Wtabs -ffixed-form -fno-align-commons
LDFLAGS := 

# Build modes
ifdef DEBUG
    FFLAGS := $(FFLAGS_DEBUG)
endif

# Default target
.PHONY: all
all: $(TARGET)

# Create directories
$(BUILDDIR):
	@mkdir -p $(BUILDDIR)

# Build executable
$(TARGET): $(OBJECTS)
	$(FC) $(FFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)
	@echo "Built $@"


# Install target
.PHONY: install
install: $(TARGET)
	@echo "Executable $(TARGET) ready in current directory"

# Test target
.PHONY: test
test: $(TARGET)
	@echo "Running tests..."
	cd $(TESTDIR) && ../$(TARGET) < cm2lab.in

# Debug build
.PHONY: debug
debug:
	$(MAKE) DEBUG=1

# Clean targets
.PHONY: clean
clean:
	rm -f $(BUILDDIR)/*.o $(BUILDDIR)/*.mod
	rm -f *.mod core
	rm -f $(TARGET)

.PHONY: distclean
distclean: clean
	rm -rf $(BUILDDIR)
	rm -f $(TARGET)

# Show help
.PHONY: help
help:
	@echo "Available targets:"
	@echo "  all       - Build the executable (default)"
	@echo "  debug     - Build with debug flags"
	@echo "  test      - Run tests"
	@echo "  install   - Install executable"
	@echo "  clean     - Remove object files and modules"
	@echo "  distclean - Remove all generated files"
	@echo "  help      - Show this help"
	@echo ""
	@echo "Variables:"
	@echo "  FC        - Fortran compiler (default: gfortran)"
	@echo "  DEBUG=1   - Enable debug build"

# Dependency tracking - explicit dependencies for each object file
$(BUILDDIR)/cm2lab.o: cm2lab.f | $(BUILDDIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILDDIR)/gauss.o: gauss.f90 | $(BUILDDIR)
	$(FC) $(FFLAGS) -c $< -o $@