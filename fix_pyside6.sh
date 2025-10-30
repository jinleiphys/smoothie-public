#!/bin/bash

# SMOOTHIE GUI - PySide6 Installation Fix Script
# Use this if you get "ModuleNotFoundError: No module named 'PySide6'"

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo "=================================================="
echo "  SMOOTHIE GUI - PySide6 Fix Script"
echo "=================================================="
echo ""

# Get conda base path
CONDA_BASE=$(conda info --base 2>/dev/null)
if [ -z "$CONDA_BASE" ]; then
    echo -e "${RED}[ERROR]${NC} conda not found!"
    exit 1
fi

# Activate environment
ENV_NAME="smoothie_gui"
echo -e "${BLUE}[INFO]${NC} Activating conda environment: $ENV_NAME"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate $ENV_NAME

echo ""
echo -e "${BLUE}[INFO]${NC} Diagnostic Information:"
echo "  Conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"
echo "  Python executable: $(which python)"
echo "  Python version: $(python --version)"
echo "  pip executable: $(which pip)"
echo ""

# Check if PySide6 is already installed
echo -e "${BLUE}[INFO]${NC} Checking current PySide6 installation..."
if python -c "import PySide6.QtWidgets; import PySide6; print('PySide6 version:', PySide6.__version__)" 2>/dev/null; then
    echo -e "${GREEN}[SUCCESS]${NC} PySide6 is already installed and working!"
    echo ""
    echo "Try running the GUI again:"
    echo "  cd smoothie_gui"
    echo "  python main.py"
    exit 0
fi

echo -e "${YELLOW}[WARNING]${NC} PySide6 not found or not working"
echo ""

# Method 1: Upgrade pip and install via pip
echo -e "${BLUE}[INFO]${NC} Method 1: Installing PySide6 via pip..."
pip install --upgrade pip
pip install --upgrade --force-reinstall PySide6>=6.5.0

if python -c "import PySide6.QtWidgets" 2>/dev/null; then
    echo -e "${GREEN}[SUCCESS]${NC} PySide6 installed successfully via pip!"
    python -c "import PySide6; print('PySide6 version:', PySide6.__version__)"
    echo ""
    echo "You can now run the GUI:"
    echo "  cd smoothie_gui"
    echo "  python main.py"
    exit 0
fi

# Method 2: Try conda installation
echo -e "${YELLOW}[WARNING]${NC} pip installation failed, trying conda..."
echo -e "${BLUE}[INFO]${NC} Method 2: Installing PySide6 via conda..."
conda install -y -c conda-forge pyside6

if python -c "import PySide6.QtWidgets" 2>/dev/null; then
    echo -e "${GREEN}[SUCCESS]${NC} PySide6 installed successfully via conda!"
    python -c "import PySide6; print('PySide6 version:', PySide6.__version__)"
    echo ""
    echo "You can now run the GUI:"
    echo "  cd smoothie_gui"
    echo "  python main.py"
    exit 0
fi

# If we get here, both methods failed
echo ""
echo -e "${RED}[ERROR]${NC} Failed to install PySide6 with both methods!"
echo ""
echo "Please try the following manually:"
echo ""
echo "  1. Check your conda environment:"
echo "     conda activate smoothie_gui"
echo "     which python"
echo "     which pip"
echo ""
echo "  2. Try installing dependencies manually:"
echo "     pip install --upgrade pip setuptools wheel"
echo "     pip install PySide6"
echo ""
echo "  3. If on macOS, you may need Xcode Command Line Tools:"
echo "     xcode-select --install"
echo ""
echo "  4. Check Python packages:"
echo "     pip list | grep -i pyside"
echo ""
echo "  5. Try recreating the environment:"
echo "     conda deactivate"
echo "     conda env remove -n smoothie_gui"
echo "     ./setup_gui.sh"
echo ""

exit 1
