#!/bin/bash

# Auto Commit Script for SMOOTHIE
# Compiles the code and commits changes

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Function to print colored messages
print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_step() {
    echo -e "${CYAN}==>${NC} $1"
}

# Get the directory where the script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

echo "=================================================="
echo "  SMOOTHIE Auto Commit Script"
echo "=================================================="
echo ""

# Step 1: Detect and select Fortran compiler
print_step "Checking for Fortran compiler..."
echo ""

IFORT_FOUND=false
GFORTRAN_FOUND=false
SELECTED_COMPILER=""

# Check for ifort
if command -v ifort &> /dev/null; then
    IFORT_FOUND=true
    IFORT_VERSION=$(ifort --version 2>&1 | head -n 1)
    print_success "Intel Fortran (ifort) found: $IFORT_VERSION"
fi

# Check for gfortran
if command -v gfortran &> /dev/null; then
    GFORTRAN_FOUND=true
    GFORTRAN_VERSION=$(gfortran --version | head -n 1)
    print_success "GNU Fortran (gfortran) found: $GFORTRAN_VERSION"
fi

# Compiler selection logic
if [ "$IFORT_FOUND" = true ] && [ "$GFORTRAN_FOUND" = true ]; then
    # Both compilers found - ask user to choose
    echo ""
    print_warning "Both Intel Fortran and GNU Fortran are available!"
    echo ""
    echo "  1) Use ifort (Intel Fortran)"
    echo "  2) Use gfortran (GNU Fortran)"
    echo ""
    read -p "Select compiler [1-2]: " -n 1 -r COMPILER_CHOICE
    echo ""

    case $COMPILER_CHOICE in
        1)
            SELECTED_COMPILER="ifort"
            print_info "Using Intel Fortran (ifort)"
            ;;
        2)
            SELECTED_COMPILER="gfortran"
            print_info "Using GNU Fortran (gfortran)"
            ;;
        *)
            print_error "Invalid choice. Defaulting to gfortran."
            SELECTED_COMPILER="gfortran"
            ;;
    esac

elif [ "$IFORT_FOUND" = true ]; then
    # Only ifort found
    SELECTED_COMPILER="ifort"
    print_info "Using Intel Fortran (ifort)"

elif [ "$GFORTRAN_FOUND" = true ]; then
    # Only gfortran found
    SELECTED_COMPILER="gfortran"
    print_info "Using GNU Fortran (gfortran)"

else
    # No compiler found - need to install gfortran
    print_warning "No Fortran compiler found!"
    echo ""
    read -p "Do you want to install gfortran? (y/N): " -n 1 -r
    echo ""

    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_error "Fortran compiler is required to compile SMOOTHIE. Exiting."
        exit 1
    fi

    print_info "Installing gfortran..."

    # Check for conda first (preferred method)
    if command -v conda &> /dev/null; then
        print_info "Installing gfortran via conda..."

        # Get conda base path and source it
        CONDA_BASE=$(conda info --base 2>/dev/null)
        if [ -n "$CONDA_BASE" ]; then
            source "$CONDA_BASE/etc/profile.d/conda.sh"
        fi

        conda install -y -c conda-forge gfortran
        print_success "gfortran installed via conda successfully!"

    # Check for package managers
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        # macOS
        if command -v brew &> /dev/null; then
            print_info "Installing gfortran via Homebrew..."
            brew install gcc
            print_success "gfortran installed via Homebrew successfully!"
        else
            print_error "Homebrew not found. Please install Homebrew first:"
            echo "  /bin/bash -c \"\$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)\""
            exit 1
        fi

    elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
        # Linux
        if command -v apt-get &> /dev/null; then
            print_info "Installing gfortran via apt-get..."
            sudo apt-get update && sudo apt-get install -y gfortran
            print_success "gfortran installed successfully!"
        elif command -v yum &> /dev/null; then
            print_info "Installing gfortran via yum..."
            sudo yum install -y gcc-gfortran
            print_success "gfortran installed successfully!"
        else
            print_error "No package manager found. Please install gfortran manually."
            exit 1
        fi

    else
        print_error "Unsupported operating system: $OSTYPE"
        exit 1
    fi

    # Verify installation
    if command -v gfortran &> /dev/null; then
        SELECTED_COMPILER="gfortran"
        GFORTRAN_VERSION=$(gfortran --version | head -n 1)
        print_success "gfortran is now available: $GFORTRAN_VERSION"
    else
        print_error "gfortran installation failed. Please install manually."
        exit 1
    fi
fi

# Update make.inc with selected compiler
if [ -n "$SELECTED_COMPILER" ]; then
    print_info "Configuring make.inc for $SELECTED_COMPILER..."

    if [ ! -f "make.inc" ]; then
        if [ "$SELECTED_COMPILER" = "ifort" ] && [ -f "make.inc.ifort" ]; then
            cp make.inc.ifort make.inc
            print_success "Created make.inc from make.inc.ifort"
        elif [ -f "make.inc.gfortran" ]; then
            cp make.inc.gfortran make.inc
            print_success "Created make.inc from make.inc.gfortran"
        else
            print_error "Template make.inc file not found!"
            exit 1
        fi
    else
        # Update existing make.inc
        if [ "$SELECTED_COMPILER" = "ifort" ]; then
            sed -i.bak 's/^FC = .*/FC = ifort/' make.inc
            sed -i.bak 's/^F90 = .*/F90 = ifort/' make.inc
        else
            sed -i.bak 's/^FC = .*/FC = gfortran/' make.inc
            sed -i.bak 's/^F90 = .*/F90 = gfortran/' make.inc
        fi
        print_success "Updated make.inc to use $SELECTED_COMPILER"
    fi
fi

echo ""

# Step 2: Check and configure LAPACK libraries
print_step "Checking for LAPACK libraries..."
LAPACK_FOUND=false
LAPACK_SOURCE=""

# Check for LAPACK in common locations
if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS: Check multiple sources in order of preference

    # 1. Check Homebrew installation
    if [ -d "/opt/homebrew/opt/lapack" ] && [ -f "/opt/homebrew/opt/lapack/lib/liblapack.dylib" ]; then
        LAPACK_FOUND=true
        LAPACK_SOURCE="Homebrew (/opt/homebrew/opt/lapack)"
    elif [ -d "/usr/local/opt/lapack" ] && [ -f "/usr/local/opt/lapack/lib/liblapack.dylib" ]; then
        LAPACK_FOUND=true
        LAPACK_SOURCE="Homebrew (/usr/local/opt/lapack)"
    # 2. Check MacPorts installation
    elif [ -f "/opt/local/lib/liblapack.dylib" ]; then
        LAPACK_FOUND=true
        LAPACK_SOURCE="MacPorts (/opt/local/lib)"
    # 3. Check Accelerate framework
    elif [ -f "/System/Library/Frameworks/Accelerate.framework/Accelerate" ] || \
         [ -f "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework/Accelerate" ]; then
        LAPACK_FOUND=true
        LAPACK_SOURCE="macOS Accelerate framework"
    fi

    if [ "$LAPACK_FOUND" = true ]; then
        print_success "LAPACK found (via $LAPACK_SOURCE)"
    fi
else
    # Linux: check for liblapack
    if ldconfig -p 2>/dev/null | grep -q liblapack; then
        LAPACK_FOUND=true
        LAPACK_SOURCE="system ldconfig"
        print_success "LAPACK found (via $LAPACK_SOURCE)"
    elif [ -f "/usr/lib/liblapack.so" ] || \
         [ -f "/usr/lib64/liblapack.so" ] || \
         [ -f "/usr/lib/x86_64-linux-gnu/liblapack.so" ]; then
        LAPACK_FOUND=true
        LAPACK_SOURCE="system libraries"
        print_success "LAPACK found (via $LAPACK_SOURCE)"
    fi
fi

# If not found, offer to install via conda
if [ "$LAPACK_FOUND" = false ]; then
    print_warning "LAPACK libraries not found!"
    echo ""
    read -p "Do you want to install LAPACK via conda? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        print_info "Installing LAPACK..."

        # Initialize conda if not already done
        if ! command -v conda &> /dev/null; then
            source "$HOME/miniconda3/etc/profile.d/conda.sh"
        fi

        # Install LAPACK and BLAS via conda
        print_info "Installing LAPACK/BLAS via conda-forge..."
        conda install -y -c conda-forge openblas lapack

        print_success "LAPACK installed successfully!"
        LAPACK_FOUND=true
        LAPACK_SOURCE="conda-forge"
    else
        print_warning "LAPACK is required for optimal performance. Continuing anyway..."
        print_info "You may need to install LAPACK manually later."
    fi
fi
echo ""

# Configure make.inc with LAPACK paths
print_info "Configuring make.inc for your system..."

if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS: Configure based on detected LAPACK source
    # Check if Accelerate framework is available
    ACCELERATE_AVAILABLE=false
    if [ -d "/System/Library/Frameworks/Accelerate.framework" ]; then
        ACCELERATE_AVAILABLE=true
    fi

    if [[ "$LAPACK_SOURCE" == *"Homebrew"* ]]; then
        # Homebrew LAPACK detected
        if [ -d "/opt/homebrew/opt/lapack" ]; then
            LAPACK_PATH="/opt/homebrew/opt/lapack/lib"
        else
            LAPACK_PATH="/usr/local/opt/lapack/lib"
        fi

        # On macOS, Accelerate framework is highly optimized by Apple
        if [ "$ACCELERATE_AVAILABLE" = true ]; then
            print_info "Both Homebrew LAPACK and Accelerate framework are available"
            print_info "Using Accelerate framework (recommended for macOS)"
            LAPACK_LIB="-framework Accelerate"
        else
            LAPACK_LIB="-L$LAPACK_PATH -llapack -lblas"
        fi
    elif [[ "$LAPACK_SOURCE" == *"MacPorts"* ]]; then
        # MacPorts LAPACK detected
        if [ "$ACCELERATE_AVAILABLE" = true ]; then
            print_info "Both MacPorts LAPACK and Accelerate framework are available"
            print_info "Using Accelerate framework (recommended for macOS)"
            LAPACK_LIB="-framework Accelerate"
        else
            LAPACK_LIB="-L/opt/local/lib -llapack"
        fi
    elif [ "$ACCELERATE_AVAILABLE" = true ]; then
        # Use Accelerate framework only
        LAPACK_LIB="-framework Accelerate"
        print_info "Using macOS Accelerate framework"
    else
        print_warning "No LAPACK library found on macOS"
        LAPACK_LIB=""
    fi
else
    # Linux: Check for conda-installed LAPACK first
    CONDA_BASE=$(conda info --base 2>/dev/null)
    if [ -n "$CONDA_BASE" ] && [ -f "$CONDA_BASE/lib/libopenblas.so" ]; then
        LAPACK_LIB="-L$CONDA_BASE/lib -lopenblas -llapack"
    elif [ -f "/usr/lib/liblapack.so" ] || [ -f "/usr/lib64/liblapack.so" ]; then
        LAPACK_LIB="-llapack -lblas"
    elif [ -f "/usr/lib/x86_64-linux-gnu/liblapack.so" ]; then
        LAPACK_LIB="-llapack -lblas"
    else
        LAPACK_LIB="-llapack"
    fi
fi

# Update LIBSTD1 in make.inc
sed -i.bak "s|^LIBSTD1 =.*|LIBSTD1 = $LAPACK_LIB|" make.inc
print_success "Updated make.inc with LAPACK configuration: $LAPACK_LIB"

echo ""

# Step 3: Select compilation mode
print_step "Select compilation mode:"
echo ""
echo "  1) Compile without GUI"
echo "  2) Compile with GUI setup"
echo ""
read -p "Enter your choice [1-2]: " -n 1 -r COMPILE_CHOICE
echo ""
echo ""

# Step 4: Compile based on choice
if [ "$COMPILE_CHOICE" = "2" ]; then
    print_info "Running GUI setup script..."
    bash setup_gui.sh
    print_success "GUI setup completed!"
    exit 0
fi

# Compile without GUI
print_step "Compiling SMOOTHIE and cm2lab..."
echo ""

# Build general modules
print_info "Building general modules..."
cd general_modules
make clean 2>/dev/null || true
make

# Build mesh modules
print_info "Building mesh modules..."
cd ../mesh_modules
make clean 2>/dev/null || true
make

# Build pw modules
print_info "Building partial wave modules..."
cd ../pw_modules
make clean 2>/dev/null || true
make

# Build pot modules
print_info "Building potential modules..."
cd ../pot_modules
make clean 2>/dev/null || true
make

# Build smoothie
print_info "Building SMOOTHIE main program..."
cd ../smoothie
make clean 2>/dev/null || true
make

# Build cm2lab
print_info "Building cm2lab utility..."
cd ../cm2lab
make clean 2>/dev/null || true
make

cd "$SCRIPT_DIR"

echo ""
echo "=================================================="
print_success "Compilation Completed Successfully!"
echo "=================================================="
echo ""
print_info "Executables are located at:"
echo ""
echo -e "  ${GREEN}✓${NC} smoothie/smoothie"
echo -e "  ${GREEN}✓${NC} cm2lab/cm2lab"
echo ""
print_info "To run SMOOTHIE:"
echo -e "  ${CYAN}cd smoothie && ./smoothie < input.in${NC}"
echo ""
print_info "To run CM2LAB coordinate conversion:"
echo -e "  ${CYAN}cd cm2lab && ./cm2lab < cm2lab.in${NC}"
echo ""
