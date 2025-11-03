#!/bin/bash

# SMOOTHIE GUI Setup Script
# This script sets up the conda environment, installs dependencies, and builds SMOOTHIE

set -e  # Exit on error

echo "=================================================="
echo "  SMOOTHIE GUI Setup Script"
echo "=================================================="
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
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

# Get the directory where the script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

print_info "Working directory: $SCRIPT_DIR"
echo ""

# Step 1: Check and install conda if needed
print_info "Checking for conda installation..."
if ! command -v conda &> /dev/null; then
    print_warning "conda not found!"
    echo ""
    read -p "Do you want to install Miniconda? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        print_info "Installing Miniconda..."

        # Detect OS
        if [[ "$OSTYPE" == "darwin"* ]]; then
            # macOS
            if [[ $(uname -m) == "arm64" ]]; then
                MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh"
            else
                MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
            fi
        elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
            # Linux
            MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
        else
            print_error "Unsupported operating system: $OSTYPE"
            exit 1
        fi

        # Download and install Miniconda
        TEMP_INSTALLER="/tmp/miniconda_installer.sh"
        print_info "Downloading Miniconda from $MINICONDA_URL..."

        # Check for download tools and use the first available one
        if command -v curl &> /dev/null; then
            curl -L -o "$TEMP_INSTALLER" "$MINICONDA_URL"
        elif command -v wget &> /dev/null; then
            wget -O "$TEMP_INSTALLER" "$MINICONDA_URL"
        else
            print_error "Neither curl nor wget found!"
            print_error "Please install one of them first:"
            echo "  Ubuntu/Debian: sudo apt-get install curl"
            echo "  CentOS/RHEL:   sudo yum install curl"
            echo "  macOS:         brew install curl"
            exit 1
        fi

        print_info "Running Miniconda installer..."
        bash "$TEMP_INSTALLER" -b -p "$HOME/miniconda3"
        rm "$TEMP_INSTALLER"

        # Initialize conda
        print_info "Initializing conda..."
        "$HOME/miniconda3/bin/conda" init bash

        # Source conda for current session
        source "$HOME/miniconda3/etc/profile.d/conda.sh"

        # Accept conda Terms of Service automatically
        print_info "Accepting conda Terms of Service..."
        conda config --set channel_priority flexible 2>/dev/null || true
        conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main 2>/dev/null || true
        conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r 2>/dev/null || true

        print_success "Miniconda installed successfully!"
        print_warning "Please restart your terminal or run: source ~/.bashrc (or ~/.zshrc)"
    else
        print_error "conda is required to continue. Exiting."
        exit 1
    fi
else
    print_success "conda found: $(conda --version)"

    # Accept conda Terms of Service (in case not already accepted)
    print_info "Ensuring conda Terms of Service are accepted..."
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main 2>/dev/null || true
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r 2>/dev/null || true
fi
echo ""

# Step 2: Check and install gfortran if needed
print_info "Checking for Fortran compiler..."
if ! command -v gfortran &> /dev/null; then
    print_warning "gfortran not found!"
    echo ""
    read -p "Do you want to install gfortran? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        print_info "Installing gfortran..."

        # Initialize conda if not already done
        if ! command -v conda &> /dev/null; then
            source "$HOME/miniconda3/etc/profile.d/conda.sh"
        fi

        # Install gfortran via conda
        print_info "Installing gfortran via conda-forge..."
        conda install -y -c conda-forge gfortran

        print_success "gfortran installed successfully!"
    else
        print_error "gfortran is required to compile SMOOTHIE. Exiting."
        exit 1
    fi
else
    print_success "Fortran compiler found: $(gfortran --version | head -n 1)"
fi
echo ""

# Step 3: Check and install LAPACK if needed
print_info "Checking for LAPACK libraries..."
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
    else
        print_warning "LAPACK is required for optimal performance. Continuing anyway..."
        print_info "You may need to install LAPACK manually later."
    fi
fi
echo ""

# Step 4: Create or update conda environment
ENV_NAME="smoothie_gui"
print_info "Setting up conda environment: $ENV_NAME"

if conda env list | grep -q "^$ENV_NAME "; then
    print_warning "Environment $ENV_NAME already exists"
    read -p "Do you want to recreate it? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        print_info "Removing existing environment..."
        conda env remove -n $ENV_NAME -y
        print_info "Creating new environment..."
        conda create -n $ENV_NAME python=3.10 -y
    else
        print_info "Using existing environment"
    fi
else
    print_info "Creating new environment..."
    conda create -n $ENV_NAME python=3.10 -y
fi

print_success "Conda environment ready"
echo ""

# Step 5: Activate environment and install Python packages
print_info "Installing Python dependencies..."

# Get conda base path
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate $ENV_NAME

# Verify we're using the correct Python and pip
print_info "Using Python: $(which python)"
print_info "Using pip: $(which pip)"
print_info "Python version: $(python --version)"

# Check if pip is from conda environment (critical check!)
PIP_PATH=$(which pip)

# Smart detection: Check if pip is from the active conda environment
if [ -n "$CONDA_PREFIX" ]; then
    # If pip path contains the conda environment prefix, it's the correct pip
    if [[ "$PIP_PATH" == "$CONDA_PREFIX"* ]]; then
        print_success "pip is correctly from conda environment"
        PIP_CMD="pip"
    else
        print_warning "pip is NOT from the conda environment!"
        print_warning "  Conda environment: $CONDA_PREFIX"
        print_warning "  pip location: $PIP_PATH"
        print_info "Using 'python -m pip' instead to ensure correct environment"
        PIP_CMD="python -m pip"
    fi
else
    print_warning "CONDA_PREFIX not set, using 'python -m pip' for safety"
    PIP_CMD="python -m pip"
fi

# Install requirements
cd smoothie_gui

# Upgrade pip first to avoid installation issues
print_info "Upgrading pip..."
$PIP_CMD install --upgrade pip

print_info "Installing requirements from requirements.txt..."
$PIP_CMD install -r requirements.txt

# Verify PySide6 installation
print_info "Verifying PySide6 installation..."
if python -c "import PySide6.QtWidgets; print('PySide6 version:', PySide6.__version__)" 2>/dev/null; then
    print_success "PySide6 installed and verified successfully"
else
    print_warning "PySide6 installation verification failed!"
    print_info "Attempting to reinstall PySide6 with verbose output..."

    # Try installing with verbose output to see what's wrong
    $PIP_CMD install --upgrade --force-reinstall --verbose PySide6>=6.5.0

    # Try again
    if python -c "import PySide6.QtWidgets; print('PySide6 version:', PySide6.__version__)" 2>/dev/null; then
        print_success "PySide6 reinstalled successfully"
    else
        print_error "Failed to install PySide6!"
        print_info "Diagnostic information:"
        echo "  Python executable: $(which python)"
        echo "  pip executable: $(which pip)"
        echo "  Python site-packages:"
        python -c "import site; print('  ', site.getsitepackages())"
        echo ""
        print_info "Trying alternative installation methods..."

        # Try conda installation as fallback
        print_info "Attempting conda installation of PySide6..."
        conda install -y -c conda-forge pyside6

        if python -c "import PySide6.QtWidgets" 2>/dev/null; then
            print_success "PySide6 installed via conda successfully"
        else
            print_error "All installation methods failed. Please try manually:"
            echo "  conda activate smoothie_gui"
            echo "  pip install --upgrade pip"
            echo "  pip install PySide6"
            echo "  # OR try: conda install -c conda-forge pyside6"
            exit 1
        fi
    fi
fi

print_success "Python dependencies installed"
echo ""

# Step 6: Build SMOOTHIE
print_info "Building SMOOTHIE Fortran code..."
cd "$SCRIPT_DIR"

# Check if make.inc exists, if not create it
if [ ! -f "make.inc" ]; then
    print_info "Creating make.inc from template..."
    if [ -f "make.inc.gfortran" ]; then
        cp make.inc.gfortran make.inc
        print_success "Created make.inc from make.inc.gfortran"
    else
        print_error "make.inc.gfortran not found!"
        exit 1
    fi
fi

# Update make.inc with correct LAPACK paths
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
        # Give user choice between Homebrew LAPACK and Accelerate
        if [ "$ACCELERATE_AVAILABLE" = true ]; then
            print_info "Both Homebrew LAPACK and Accelerate framework are available"
            print_info "Accelerate framework is Apple's optimized implementation"
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
print_success "SMOOTHIE compiled successfully"
echo ""

# Step 7: Create launcher script
print_info "Creating launcher script..."
cat > run_smoothie_gui.sh << 'EOF'
#!/bin/bash

# SMOOTHIE GUI Launcher Script

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Get conda base path
CONDA_BASE=$(conda info --base 2>/dev/null)
if [ -z "$CONDA_BASE" ]; then
    echo "Error: conda not found!"
    exit 1
fi

# Activate environment
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate smoothie_gui

# Launch GUI
cd "$SCRIPT_DIR/smoothie_gui"
python main.py
EOF

chmod +x run_smoothie_gui.sh
print_success "Launcher script created: run_smoothie_gui.sh"
echo ""

# Step 8: Print completion message
echo "=================================================="
print_success "Setup completed successfully!"
echo "=================================================="
echo ""
echo "To run SMOOTHIE GUI:"
echo ""
echo -e "  ${BLUE}Recommended - Use the launcher script:${NC}"
echo -e "    ${GREEN}./run_smoothie_gui.sh${NC}"
echo ""
echo -e "  ${BLUE}Alternative - Manual launch:${NC}"
echo -e "    ${GREEN}conda activate $ENV_NAME${NC}"
echo -e "    ${GREEN}cd smoothie_gui${NC}"
echo -e "    ${GREEN}python main.py${NC}"
echo ""
print_warning "IMPORTANT: You MUST activate the conda environment before running!"
echo -e "  ${RED}DON'T run:${NC} python main.py  (without activating environment)"
echo -e "  ${GREEN}DO run:${NC}    conda activate smoothie_gui && cd smoothie_gui && python main.py"
echo ""
print_info "The GUI provides:"
echo "  - Parameter input forms with validation"
echo "  - Real-time output monitoring"
echo "  - Integrated result plotting"
echo "  - Light and dark themes"
echo "  - File save/load functionality"
echo ""
print_info "For help and documentation, see smoothie_gui/README.md"
echo ""
