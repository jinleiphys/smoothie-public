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

# Step 1: Check if conda is installed
print_info "Checking for conda installation..."
if ! command -v conda &> /dev/null; then
    print_error "conda not found! Please install Anaconda or Miniconda first."
    echo "Download from: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi
print_success "conda found: $(conda --version)"
echo ""

# Step 2: Check if gfortran is installed
print_info "Checking for Fortran compiler..."
if ! command -v gfortran &> /dev/null; then
    print_warning "gfortran not found!"
    print_info "Installing gfortran via brew (macOS)..."
    if command -v brew &> /dev/null; then
        brew install gcc
    else
        print_error "Please install gfortran manually or install Homebrew first."
        exit 1
    fi
fi
print_success "Fortran compiler found: $(gfortran --version | head -n 1)"
echo ""

# Step 3: Create or update conda environment
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

# Step 4: Activate environment and install Python packages
print_info "Installing Python dependencies..."

# Get conda base path
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate $ENV_NAME

# Install requirements
cd smoothie_gui
pip install -r requirements.txt

print_success "Python dependencies installed"
echo ""

# Step 5: Build SMOOTHIE
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

# Step 6: Create launcher script
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

# Step 7: Print completion message
echo "=================================================="
print_success "Setup completed successfully!"
echo "=================================================="
echo ""
echo "To run SMOOTHIE GUI:"
echo ""
echo "  Option 1 (Quick launch):"
echo -e "    ${GREEN}./run_smoothie_gui.sh${NC}"
echo ""
echo "  Option 2 (Manual):"
echo -e "    ${GREEN}conda activate $ENV_NAME${NC}"
echo -e "    ${GREEN}cd smoothie_gui${NC}"
echo -e "    ${GREEN}python main.py${NC}"
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
