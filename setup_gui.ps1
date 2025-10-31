# SMOOTHIE GUI Setup Script for Windows
# This script sets up the conda environment, installs dependencies, and builds SMOOTHIE

$ErrorActionPreference = "Stop"

Write-Host "==================================================" -ForegroundColor Cyan
Write-Host "  SMOOTHIE GUI Setup Script (Windows)" -ForegroundColor Cyan
Write-Host "==================================================" -ForegroundColor Cyan
Write-Host ""

# Functions for colored output
function Write-Info {
    param($Message)
    Write-Host "[INFO] $Message" -ForegroundColor Blue
}

function Write-Success {
    param($Message)
    Write-Host "[SUCCESS] $Message" -ForegroundColor Green
}

function Write-Warning {
    param($Message)
    Write-Host "[WARNING] $Message" -ForegroundColor Yellow
}

function Write-Error {
    param($Message)
    Write-Host "[ERROR] $Message" -ForegroundColor Red
}

# Get script directory
$SCRIPT_DIR = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $SCRIPT_DIR

Write-Info "Working directory: $SCRIPT_DIR"
Write-Host ""

# Step 1: Check and install conda if needed
Write-Info "Checking for conda installation..."
$condaExists = $null -ne (Get-Command conda -ErrorAction SilentlyContinue)

if (-not $condaExists) {
    Write-Warning "conda not found!"
    Write-Host ""
    $install = Read-Host "Do you want to install Miniconda? (y/N)"

    if ($install -match "^[Yy]$") {
        Write-Info "Installing Miniconda..."

        # Download Miniconda for Windows
        $MINICONDA_URL = "https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe"
        $TEMP_INSTALLER = "$env:TEMP\miniconda_installer.exe"

        Write-Info "Downloading Miniconda from $MINICONDA_URL..."
        Invoke-WebRequest -Uri $MINICONDA_URL -OutFile $TEMP_INSTALLER

        Write-Info "Running Miniconda installer..."
        Write-Host "Please follow the installer prompts. It's recommended to:"
        Write-Host "  - Install for 'Just Me'"
        Write-Host "  - Add Miniconda to PATH (check the box)"
        Write-Host "  - Register Miniconda as default Python"
        Start-Process -FilePath $TEMP_INSTALLER -Wait

        Remove-Item $TEMP_INSTALLER -ErrorAction SilentlyContinue

        Write-Success "Miniconda installed successfully!"
        Write-Warning "Please restart PowerShell and run this script again."
        Write-Host ""
        Read-Host "Press Enter to exit"
        exit 0
    } else {
        Write-Error "conda is required to continue. Exiting."
        Read-Host "Press Enter to exit"
        exit 1
    }
} else {
    $condaVersion = (conda --version 2>&1)
    Write-Success "conda found: $condaVersion"
}
Write-Host ""

# Initialize conda for PowerShell
Write-Info "Initializing conda for PowerShell..."
try {
    $condaPath = (Get-Command conda).Source
    $condaRoot = Split-Path (Split-Path $condaPath)
    $condaHook = Join-Path $condaRoot "shell\condabin\conda-hook.ps1"

    if (Test-Path $condaHook) {
        . $condaHook
    }
} catch {
    Write-Warning "Could not initialize conda hook. Continuing anyway..."
}

# Step 2: Check and install gfortran if needed
Write-Info "Checking for Fortran compiler..."
$gfortranExists = $null -ne (Get-Command gfortran -ErrorAction SilentlyContinue)

if (-not $gfortranExists) {
    Write-Warning "gfortran not found!"
    Write-Host ""
    $install = Read-Host "Do you want to install gfortran via conda? (y/N)"

    if ($install -match "^[Yy]$") {
        Write-Info "Installing gfortran..."
        conda install -y -c conda-forge m2w64-gcc-fortran
        Write-Success "gfortran installed successfully!"
    } else {
        Write-Warning "gfortran is required to compile SMOOTHIE."
        Write-Info "You can also install it via MinGW-w64 or WSL."
        Write-Host ""
        $continue = Read-Host "Continue anyway? (y/N)"
        if ($continue -notmatch "^[Yy]$") {
            exit 1
        }
    }
} else {
    $gfortranVersion = (gfortran --version 2>&1 | Select-Object -First 1)
    Write-Success "Fortran compiler found: $gfortranVersion"
}
Write-Host ""

# Step 3: Create or update conda environment
$ENV_NAME = "smoothie_gui"
Write-Info "Setting up conda environment: $ENV_NAME"

$envExists = (conda env list | Select-String -Pattern "^$ENV_NAME\s")

if ($envExists) {
    Write-Warning "Environment $ENV_NAME already exists"
    $recreate = Read-Host "Do you want to recreate it? (y/N)"

    if ($recreate -match "^[Yy]$") {
        Write-Info "Removing existing environment..."
        conda env remove -n $ENV_NAME -y
        Write-Info "Creating new environment..."
        conda create -n $ENV_NAME python=3.10 -y
    } else {
        Write-Info "Using existing environment"
    }
} else {
    Write-Info "Creating new environment..."
    conda create -n $ENV_NAME python=3.10 -y
}

Write-Success "Conda environment ready"
Write-Host ""

# Step 4: Activate environment and install Python packages
Write-Info "Installing Python dependencies..."
conda activate $ENV_NAME

Write-Info "Python version: $(python --version)"
Write-Info "Python location: $(Get-Command python).Source"

# Install requirements
Set-Location "$SCRIPT_DIR\smoothie_gui"

Write-Info "Upgrading pip..."
python -m pip install --upgrade pip

Write-Info "Installing requirements from requirements.txt..."
python -m pip install -r requirements.txt

# Verify PySide6 installation
Write-Info "Verifying PySide6 installation..."
$pyside6Test = python -c "import PySide6.QtWidgets; print('PySide6 version:', PySide6.__version__)" 2>&1

if ($LASTEXITCODE -eq 0) {
    Write-Success "PySide6 installed and verified successfully"
} else {
    Write-Warning "PySide6 installation verification failed!"
    Write-Info "Attempting to reinstall PySide6..."

    python -m pip install --upgrade --force-reinstall PySide6

    $pyside6Test2 = python -c "import PySide6.QtWidgets" 2>&1
    if ($LASTEXITCODE -eq 0) {
        Write-Success "PySide6 reinstalled successfully"
    } else {
        Write-Error "Failed to install PySide6!"
        Write-Info "Trying conda installation as fallback..."
        conda install -y -c conda-forge pyside6

        $pyside6Test3 = python -c "import PySide6.QtWidgets" 2>&1
        if ($LASTEXITCODE -ne 0) {
            Write-Error "All installation methods failed."
            Write-Host "Please try manually:"
            Write-Host "  conda activate smoothie_gui"
            Write-Host "  pip install PySide6"
            Read-Host "Press Enter to exit"
            exit 1
        }
    }
}

Write-Success "Python dependencies installed"
Write-Host ""

# Step 5: Build SMOOTHIE
Write-Info "Building SMOOTHIE Fortran code..."
Set-Location $SCRIPT_DIR

# Check if make.inc exists
if (-not (Test-Path "make.inc")) {
    Write-Info "Creating make.inc from template..."
    if (Test-Path "make.inc.gfortran") {
        Copy-Item "make.inc.gfortran" "make.inc"
        Write-Success "Created make.inc from make.inc.gfortran"
    } else {
        Write-Error "make.inc.gfortran not found!"
        Read-Host "Press Enter to exit"
        exit 1
    }
}

Write-Info "Note: Building on Windows requires make utility."
Write-Info "You can install it via:"
Write-Info "  - chocolatey: choco install make"
Write-Info "  - conda: conda install -c conda-forge make"
Write-Info "  - Or use WSL (Windows Subsystem for Linux)"
Write-Host ""

$hasMake = $null -ne (Get-Command make -ErrorAction SilentlyContinue)

if ($hasMake) {
    Write-Info "Make found, building SMOOTHIE..."

    # Build modules
    Write-Info "Building general modules..."
    Set-Location "general_modules"
    make clean 2>$null
    make

    Write-Info "Building mesh modules..."
    Set-Location "..\mesh_modules"
    make clean 2>$null
    make

    Write-Info "Building partial wave modules..."
    Set-Location "..\pw_modules"
    make clean 2>$null
    make

    Write-Info "Building potential modules..."
    Set-Location "..\pot_modules"
    make clean 2>$null
    make

    Write-Info "Building SMOOTHIE main program..."
    Set-Location "..\smoothie"
    make clean 2>$null
    make

    Write-Info "Building cm2lab utility..."
    Set-Location "..\cm2lab"
    make clean 2>$null
    make

    Set-Location $SCRIPT_DIR
    Write-Success "SMOOTHIE compiled successfully"
} else {
    Write-Warning "Make utility not found!"
    Write-Info "SMOOTHIE Fortran code needs to be compiled manually."
    Write-Info "Consider using WSL (Windows Subsystem for Linux) for easier compilation."
    Write-Host ""
    Write-Info "In WSL, you can run: ./setup_gui.sh"
}
Write-Host ""

# Step 6: Create launcher script
Write-Info "Creating launcher script..."

$launcherContent = @"
@echo off
REM SMOOTHIE GUI Launcher Script for Windows

echo Activating conda environment...
call conda activate smoothie_gui

if errorlevel 1 (
    echo Error: Failed to activate conda environment
    pause
    exit /b 1
)

echo Launching SMOOTHIE GUI...
cd /d "%~dp0\smoothie_gui"
python main.py

pause
"@

Set-Content -Path "run_smoothie_gui.bat" -Value $launcherContent
Write-Success "Launcher script created: run_smoothie_gui.bat"
Write-Host ""

# Step 7: Print completion message
Write-Host "==================================================" -ForegroundColor Cyan
Write-Success "Setup completed successfully!"
Write-Host "==================================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "To run SMOOTHIE GUI:" -ForegroundColor Cyan
Write-Host ""
Write-Host "  Recommended - Double-click the launcher:" -ForegroundColor Blue
Write-Host "    run_smoothie_gui.bat" -ForegroundColor Green
Write-Host ""
Write-Host "  Alternative - Manual launch:" -ForegroundColor Blue
Write-Host "    conda activate $ENV_NAME" -ForegroundColor Green
Write-Host "    cd smoothie_gui" -ForegroundColor Green
Write-Host "    python main.py" -ForegroundColor Green
Write-Host ""
Write-Warning "IMPORTANT: You MUST activate the conda environment before running!"
Write-Host ""
Write-Info "The GUI provides:"
Write-Host "  - Parameter input forms with validation"
Write-Host "  - Real-time output monitoring"
Write-Host "  - Integrated result plotting"
Write-Host "  - Light and dark themes"
Write-Host "  - File save/load functionality"
Write-Host ""
Write-Info "For help and documentation, see smoothie_gui\README.md"
Write-Host ""

if (-not $hasMake) {
    Write-Host "==================================================" -ForegroundColor Yellow
    Write-Warning "IMPORTANT: SMOOTHIE Fortran code not compiled!"
    Write-Host "To use the GUI, you need to compile SMOOTHIE first." -ForegroundColor Yellow
    Write-Host "Recommended: Use WSL (Windows Subsystem for Linux)" -ForegroundColor Yellow
    Write-Host "  1. Install WSL: wsl --install" -ForegroundColor Yellow
    Write-Host "  2. Open WSL terminal" -ForegroundColor Yellow
    Write-Host "  3. Run: ./setup_gui.sh" -ForegroundColor Yellow
    Write-Host "==================================================" -ForegroundColor Yellow
    Write-Host ""
}

Read-Host "Press Enter to exit"
