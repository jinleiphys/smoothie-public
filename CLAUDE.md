# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

SMOOTHIE is a Modern Fortran code for calculating non-elastic breakup reactions using the Ichimura-Austern-Vincent formalism. It handles reactions of the form `a(=b+x) + A → b + B*` where the projectile `a` breaks up into detected particle `b` and undetected particle `x`.

## Architecture

### Module Structure
- **general_modules/**: Core Fortran modules (precision, constants, systems, channels)
- **mesh_modules/**: Numerical integration and mesh handling (angular, interpolation, Gaussian quadrature)
- **pw_modules/**: Partial wave calculations (Coulomb functions, spherical harmonics, Clebsch-Gordan)
- **pot_modules/**: Optical potential models (Woods-Saxon, global potentials like KD02, CH89)
- **smoothie/**: Main calculation engine with DWBA methods
- **cm2lab/**: Center-of-mass to laboratory frame conversion utility

### Key Components
- **DWBA Methods**: Five calculation methods with different physics and performance:
  - `dwba=1`: No spins, rbx variable (most common, PLM optimized)
  - `dwba=2`: No spins, rb variable  
  - `dwba=3`: With spins, rb variable
  - `dwba=4`: Enhanced Lagrange mesh method
  - `dwba=5`: Lagrange mesh with rb variable
- **Calculation Flow**: Main program (smoothie.F) → initialization → method selection → DWBA calculation modules
- **Green's Functions**: Core physics calculations in green.f with optimized matrix operations
- **Potential Models**: Support for Woods-Saxon, Gaussian, and global optical potentials (KD02, CH89)

## Build System

### Configuration Files
The build system uses a hierarchical Makefile structure with compiler-specific configurations:
- **make.inc**: Active compiler configuration (copy from templates below)
- **make.inc.gfortran**: GNU Fortran configuration template 
- **make.inc.ifort**: Intel Fortran configuration template
- **dir.inc**: Module directory definitions
- **Individual Makefiles**: Each module has its own Makefile

### Compiler Setup
```bash
# For GNU Fortran (recommended for most users)
cp make.inc.gfortran make.inc

# For Intel Fortran (better performance with MKL)
cp make.inc.ifort make.inc
```

### Key Build Variables
- `FC/F90`: Fortran compiler (gfortran or ifort)
- `FOPT`: Compilation flags with preprocessor definitions for version info
- `LIBSTD`: External libraries (LAPACK/BLAS required, Intel MKL optional)
- `COMPILE_OPT1`: Core optimization flags (-O3, -fopenmp, -fcheck=all for debugging)

### Essential Commands
```bash
# Initial setup and build
cp make.inc.gfortran make.inc  # Choose compiler configuration
cd smoothie/
make clean
make

# Full rebuild (when changing compiler or optimization flags)
make clean -C general_modules
make clean -C mesh_modules  
make clean -C pw_modules
make clean -C pot_modules
make clean
make

# Build utilities
cd cm2lab/
make
```

## Testing

### Test Structure
- **smoothie/test/**: Main test cases with input files (*.in) and reference outputs
- **smoothie/test/1/**: Additional test case directory
- **cm2lab/test/**: Tests for frame conversion utility

### Running Tests
```bash
# Basic functionality test (deuteron breakup)
cd smoothie/test/
./smoothie < test.in

# Alternative test case
cd smoothie/test/1/
./smoothie < test.in

# Specific nuclear reaction test
cd smoothie/test/
./smoothie < 11Be64Zn.in

# Frame conversion test  
cd cm2lab/test/
./cm2lab < cm2lab.in
```

### Input File Format
Uses Fortran namelist format with sections:
- **&GLOBAL**: Global parameters (energy, angular range, DWBA method)
- **&SYSTEM**: Particle properties (masses, charges, spins, binding energy)
- **&OUTGOING**: Energy parameters for outgoing particles
- **&POTENTIAL**: Optical potential definitions for each channel

## Performance Optimizations

The codebase includes two major performance optimizations targeting the most expensive calculations:

### PLM Caching (Active)
- **Target**: Associated Legendre Polynomial calculations in nested loops
- **Implementation**: Hash-based cache in `plm_cache.f` with 10,000-entry capacity
- **Impact**: 5-10x speedup for DWBA=1 calculations with >95% hit rates
- **Monitoring**: Automatic cache statistics reported at end of calculation
- **Memory**: <0.5% overhead with intelligent collision resolution

### BLAS Matrix Operations (Active)  
- **Target**: O(nr²) R-function calculations in `iavdwbarbx.f`
- **Implementation**: ZGEMV matrix-vector operations replace scalar loops
- **Impact**: 3-5x speedup for large radial meshes (nr>500)
- **Requirements**: LAPACK/BLAS libraries (already configured in Makefiles)
- **Validation**: Produces identical numerical results to original loops

## Dependencies

### Required Libraries
- **LAPACK/BLAS**: Essential for linear algebra operations
- **Modern Fortran compiler**: F95 or later (gfortran 4.8+, ifort 14+)

### Platform Support
- Linux (Ubuntu, CentOS, RHEL)
- macOS (with Xcode command line tools) 
- Windows (with MinGW or WSL)

## Development Guidelines

### Code Conventions
- Modern Fortran standards with meaningful variable names
- 2-space indentation
- Comments for complex physics calculations
- Module-based architecture

### File Organization
- `*.f90`: Main Fortran source files
- `*.f`: Legacy Fortran 77 files (minimize new usage)
- `*.F`: Fortran files with preprocessor directives
- `*.inc`: Include files for common parameters

### Physics Validation
All code changes should preserve numerical accuracy. Compare results with reference calculations when modifying core algorithms.

## Debugging and Development

### Compilation Flags
The build system includes debugging support through compiler flags:
- **Production**: Use `make.inc` with `-O3` optimization  
- **Debug mode**: Enable `-fcheck=all -g` in `COMPILE_OPT1` for array bounds checking
- **Profiling**: Add `-pg` flag for gprof profiling support

### Key Algorithms
- **Main calculation entry**: `smoothie.F` lines 36-49 select DWBA method
- **Green's functions**: Central physics in `green.f` 
- **Potential evaluation**: Method dispatch in `pot_modules/pot.f`
- **PLM calculations**: Performance-critical code in `plm_cache.f`
- **Matrix operations**: BLAS-optimized sections in `iavdwbarbx.f`

### Performance Analysis
- **PLM cache**: Monitor hit rates (should be >95% for good performance)
- **Memory usage**: Large problems scale as O(nr²) for Green's functions  
- **Bottlenecks**: For nr<100, PLM dominates; for nr>500, matrix operations dominate
- **Compiler optimization**: Intel Fortran with MKL typically 20-30% faster than gfortran