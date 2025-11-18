# 🌟 SMOOTHIE

<div align="center">

**S**cattering **M**odel of **O**ptical **O**perator **Th**eory for **I**chimura-**A**ustern-**V**incent **E**quations

*A Modern Fortran code for non-elastic breakup calculations in inclusive breakup reactions*

[![Website](https://img.shields.io/badge/Website-smoothie.fewbody.com-blue)](https://smoothie.fewbody.com)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey)](README.md)

</div>

---

## 📖 Overview

SMOOTHIE is a computer code developed by **Jin Lei** and **Antonio M. Moro** to perform non-elastic breakup calculations in inclusive breakup reactions using the formalism of Ichimura, Austern and Vincent. The code handles reactions of the form:

```
a(=b+x) + A → b + B*, where B* = (x+A)**
```

## 🚀 Quick Start

### Option 1: GUI Setup (Recommended for New Users)

Run the automated setup script for a complete GUI environment:

```bash
# One-command setup
chmod +x setup_gui.sh
./setup_gui.sh

# Launch the GUI
./run_smoothie_gui.sh
```

The GUI provides an intuitive interface with forms, real-time plotting, and output monitoring.

### Option 2: Command-Line Installation

**Prerequisites:**
- **Modern Fortran** compiler (gfortran, ifort, etc.)
- **Make** build system
- **LAPACK/BLAS** libraries (recommended)

**Installation:**
```bash
# 1. Configure build environment
vim make.inc

# 2. Build the program
cd smoothie/
make

# 3. Test installation
cd test/
./smoothie < test.in
```

## 🖥️ Graphical User Interface

SMOOTHIE includes a modern graphical user interface that provides an intuitive way to set up calculations, monitor progress, and visualize results in real-time.

### GUI Features

- **📝 Form-Based Input**: Organized parameter entry with tabs for Global, System, Outgoing, and Potential settings
- **📊 Real-Time Plotting**: Automatic visualization of calculation results with multiple plot types
- **📟 Output Monitoring**: Color-coded log display with live SMOOTHIE output
- **🎨 Modern Design**: Beautiful, responsive interface with light and dark theme support
- **💾 File Operations**: Load/save input files, load examples, export results
- **⚡ Integrated Execution**: Run SMOOTHIE directly from the GUI with real-time feedback

### Quick GUI Setup

The automated setup script handles all dependencies and configuration:

**Linux / macOS:**
```bash
# Run the setup script
chmod +x setup_gui.sh
./setup_gui.sh
```

**Windows:**
```powershell
# Run in PowerShell (as Administrator recommended)
.\setup_gui.ps1
```

> **Note for Windows users**: For the best experience, we recommend using WSL (Windows Subsystem for Linux) to run the bash setup script. Native Windows compilation requires additional tools (make, gfortran) which can be complex to set up.

The script will automatically:
1. ✅ Install conda (Miniconda) if not present
2. ✅ Install gfortran compiler if needed
3. ✅ Detect and configure LAPACK/BLAS libraries
4. ✅ Create `smoothie_gui` conda environment
5. ✅ Install Python dependencies (PySide6, matplotlib, numpy)
6. ✅ Build SMOOTHIE Fortran code with optimized settings
7. ✅ Create a launcher script (`run_smoothie_gui.sh`)

### Running the GUI

**Linux / macOS - Recommended method:**
```bash
./run_smoothie_gui.sh
```

**Windows - Recommended method:**
```cmd
run_smoothie_gui.bat
```
(Or double-click the `run_smoothie_gui.bat` file)

**Alternative method (all platforms):**
```bash
conda activate smoothie_gui
cd smoothie_gui
python main.py
```

### Using the GUI

1. **Input Parameters**: Navigate tabs to enter calculation parameters or load an example
2. **Run Calculation**: Click "Run" button (Ctrl+R) and monitor progress in the Output Log
3. **View Results**: Results automatically appear in the Plot tab when complete
4. **Save/Load**: Use File menu or shortcuts (Ctrl+S to save, Ctrl+O to open)

### GUI Keyboard Shortcuts

| Shortcut | Action |
|----------|--------|
| `Ctrl+N` | New file |
| `Ctrl+O` | Open file |
| `Ctrl+S` | Save file |
| `Ctrl+R` | Run SMOOTHIE |
| `Ctrl+.` | Stop calculation |
| `Ctrl+Q` | Quit application |

For detailed GUI documentation, see [`smoothie_gui/README.md`](smoothie_gui/README.md).

---

## 📁 Project Structure

```
smoothie-public/
├── general_modules/       # Core Fortran modules
│   ├── precision.F90     # Precision definitions (double/quadruple precision)
│   ├── constants.F90     # Physical and mathematical constants
│   ├── systems.f         # System properties and particle data
│   └── channels.f        # Channel coupling and reaction definitions
│
├── mesh_modules/         # Numerical integration and mesh handling
│   ├── mesh.f            # Mesh setup and radial grids
│   ├── angularmesh.f     # Angular momentum coupling and grids
│   ├── interpolation.f   # Interpolation routines
│   └── derivative.f      # Numerical derivatives
│
├── pw_modules/           # Partial wave calculations
│   ├── coul90.f          # Coulomb wave functions
│   ├── coulcc.f          # Coulomb functions with complex arguments
│   ├── whittaker.f       # Whittaker functions
│   ├── spharm.f          # Spherical harmonics
│   └── clebsch.f         # Clebsch-Gordan coefficients
│
├── pot_modules/          # Optical potential models
│   ├── pot.f             # Main potential interface and dispatcher
│   ├── kd02.f            # Koning-Delaroche (2002) global potential
│   ├── ch89.f            # Chapel-Hill 89 global potential
│   ├── daehnick.f        # Daehnick potential
│   └── yyq06.f           # YYQ06 potential model
│
├── smoothie/             # Main calculation engine
│   ├── smoothie.F        # Main program entry point
│   ├── input.f           # Input file parsing
│   ├── iavdwbarbx.f      # DWBA=1 method (rbx variable, no spins)
│   ├── iavzerospin.f     # DWBA=2 method (rb variable, no spins)
│   ├── iavdwba.f         # DWBA=3 method (rb variable, with spins)
│   ├── iavdwbarbx_lagrange.f  # DWBA=4 Lagrange mesh (rbx variable)
│   ├── iavdwbarb_lagrange.f   # DWBA=5 Lagrange mesh (rb variable)
│   ├── lagrange_mesh.f   # Lagrange mesh implementation
│   ├── plm_cache.f       # PLM caching optimization (5-10x speedup)
│   ├── green.f           # Green's function calculations
│   ├── bound.f           # Bound state wave functions
│   ├── scatt.f           # Scattering wave functions
│   ├── fuspot.f          # Fusion potential calculations
│   ├── zerorange.f       # Zero-range approximation
│   ├── test/             # Test cases and examples
│   │   ├── test.in       # Basic deuteron breakup test
│   │   └── 11Be64Zn.in   # 11Be + 64Zn reaction example
│   └── Makefile          # Build configuration
│
├── cm2lab/               # Center-of-mass to lab frame converter
│   ├── cm2lab.f          # Main conversion program
│   └── test/             # Conversion test cases
│
├── smoothie_gui/         # Python-based graphical user interface
│   ├── main.py           # GUI entry point
│   ├── main_window.py    # Main window implementation
│   ├── input_panel.py    # Input parameter forms
│   ├── plot_widget.py    # Real-time plotting
│   ├── log_widget.py     # Output monitoring
│   ├── runner.py         # SMOOTHIE execution handler
│   ├── path_utils.py     # Path and file utilities
│   ├── styles.py         # GUI styling and themes
│   └── requirements.txt  # Python dependencies
│
├── make.inc              # Active compiler configuration
├── make.inc.gfortran     # GNU Fortran template
├── make.inc.ifort        # Intel Fortran template
├── dir.inc               # Module directory definitions
├── compile.sh            # Automated compilation script
├── setup_gui.sh          # GUI environment setup script (Linux/macOS)
├── setup_gui.ps1         # GUI environment setup script (Windows)
└── run_smoothie_gui.sh   # GUI launcher (generated by setup_gui.sh)
```

### Key Directories

| Directory | Purpose | Key Features |
|-----------|---------|--------------|
| **general_modules/** | Foundation modules | Precision control, physical constants, system data structures |
| **mesh_modules/** | Numerical methods | Integration grids, interpolation, quadrature schemes |
| **pw_modules/** | Angular momentum | Coulomb functions, spherical harmonics, quantum coupling |
| **pot_modules/** | Interaction potentials | Woods-Saxon, global potentials (KD02, CH89), external files |
| **smoothie/** | Main engine | Five DWBA methods, Green's functions, PLM caching |
| **cm2lab/** | Utility | Frame transformation for experimental comparison |
| **smoothie_gui/** | User interface | Modern Qt-based GUI with plotting and monitoring |

### Build System Files

| File | Purpose |
|------|---------|
| `make.inc` | Active compiler configuration (copied from templates) |
| `make.inc.gfortran` | GNU Fortran configuration template |
| `make.inc.ifort` | Intel Fortran configuration template |
| `dir.inc` | Module directory paths |
| `compile.sh` | Automated build script with library detection |

---

## 📋 Input File Structure

SMOOTHIE uses Fortran namelist format with the following main sections:

### 🔧 &GLOBAL - Global Parameters

**Basic Settings**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `hcm` | Radial step size (fm) for wave function integration | 0.05 |
| `rmax` | Maximum radius (fm) for all channels | 50.0 |
| `lmax` | Maximum orbital angular momentum (λₐ, λᵦ) | - |
| `lmin` | Minimum orbital angular momentum | 0 |
| `lxmax` | Maximum lₓ for x+A channel | lmax |
| `elab` | Laboratory energy (MeV) of projectile "a" | - |

**Angular Distribution**

| Parameter | Description | Default |
|-----------|-------------|---------|
| `thmin`, `thmax` | Min/max scattering angles (degrees) | - |
| `thinc` | Angular increment (degrees) for output | - |
| `nx` | Gaussian quadrature points for angular integration | 34 |

**DWBA Methods**

| Value | Description |
|-------|-------------|
| `dwba=1` | No intrinsic spins, rbx integration variable |
| `dwba=2` | No intrinsic spins, rb integration variable |
| `dwba=3` | With spins, rb integration variable |
| `dwba=4` | Enhanced method with Lagrange mesh |
| `dwba=5` | Lagrange mesh with rb variable |

### 🎯 &SYSTEM - System Parameters

**Particle Properties**

| Parameter | Description |
|-----------|-------------|
| `namep`, `namet`, `nameb`, `namex` | Names of projectile, target, detected, and undetected particles |
| `massp`, `masst`, `massb`, `massx` | Masses (amu) of each particle |
| `zp`, `zt`, `zb`, `zx` | Charges of each particle |
| `jp`, `jt`, `jb`, `jx` | Spins of each particle |

**Bound State Properties**

| Parameter | Description |
|-----------|-------------|
| `lbx` | Orbital angular momentum of b-x within projectile bound state |
| `sbx` | Total spin coupling (jb + jx) for detected and undetected particles |
| `nodes` | Number of nodes in the b-x bound state wave function |
| `be` | Binding energy (MeV) of the b-x system (positive value) |

### 📊 &OUTGOING - Output Energy Parameters

| Parameter | Description |
|-----------|-------------|
| `ecmbmin` | Minimum energy (MeV) of outgoing particle b |
| `ecmbmax` | Maximum energy (MeV) of outgoing particle b |
| `ecmbh` | Energy step size (MeV) for output |

### ⚛️ &POTENTIAL - Optical Potentials

**System Identification (kp1)**

| Value | Description |
|-------|-------------|
| `'a'` | Projectile + target (a+A) channel |
| `'b'` | Detected + recoil (b+B) channel |
| `'t'` | Detected + target (b+A) channel |
| `'x'` | Undetected + target (x+A) channel |
| `'p'` | Detected + undetected (b+x) bound state |

**Potential Models (ptype)**

| Value | Description |
|-------|-------------|
| `1` | Woods-Saxon potential with full parameter set |
| `2` | Gaussian potential |
| `3` | Koning-Delaroche (KD02) nucleon-nucleus global potential |
| `4` | CH89 nucleon-nucleus global potential |
| `41-44` | Read potential from external files (fort.41 to fort.44) |

## 🔨 Compilation and Installation

### Build Options

| Command | Description |
|---------|-------------|
| `make` | Standard build |
| `make clean` | Clean build files |

### Supported Platforms
- 🐧 **Linux** (Ubuntu, CentOS, RHEL)
- 🍎 **macOS** (with Xcode command line tools)
- 🪟 **Windows** (with MinGW or WSL)

### Troubleshooting
- **Missing libraries**: Install LAPACK/BLAS development packages
- **Compiler errors**: Check Fortran compiler version (F95 or later required)
- **Linking issues**: Verify library paths in `make.inc`

## 📝 Examples

### 93Nb(d,pX) - Deuteron Breakup on Niobium

This example demonstrates a deuteron breakup reaction on ⁹³Nb target, analyzing the (d,p) channel.

**Input File:**

```fortran
NAMELIST
&GLOBAL      hcm=0.05  lmax=25  elab=25.5 thmin=0. thmax=180.  printf=f dwba=1
             thinc=1   nx=34 rmax=50   nr=100  lxmax=12  /
&SYSTEM     namep='d'     massp=2.       zp=1.0    jp=0. sbx=0.
            namet='93Nb'  masst=93.0     zt=41.0   jt=0.0  be=2.224
            nameb='p'     massb=1.0078        zb=1.0    jb=0.
            namex='n'     massx=1.0087   zx=0.0    jx=0.  lbx=0    nodes=1  /
&OUTGOING   ecmbmin=2 ecmbmax=30 ecmbh=1  /
&OUTGOING /
&POTENTIAL  kp1='a' ptype=1 a1=0 a2=93 rc=1.3
            uv=77.3 av=0.77 rv=1.15
            uw=6.1  aw=0.47 rw=1.33
            wd=8.4  awd=0.77 rwd=1.37
            /
&POTENTIAL  kp1='b'  ptype=4 a1=0 a2=94
           /
&POTENTIAL  kp1='x' ptype=4 a1=0 a2=93
            /
&POTENTIAL  kp1='p' ptype=2 a1=1 a2=1 rc=1.5
            uv=72.15 av=1.484
            /
&POTENTIAL  kp1='t' ptype=4 a1=0 a2=93
           /
&POTENTIAL /
```

**Key Parameters:**
- 🎯 **Incident Energy**: 25.5 MeV deuteron beam
- 📐 **Angular Range**: 0° to 180° in 1° steps
- ⚡ **Energy Range**: 2.0 to 30.0 MeV outgoing proton energy
- 🔬 **DWBA Method**: Standard method (dwba=1)
- 🧮 **Potential Models**: Woods-Saxon (a+A), CH89 global (others), Gaussian (bound state)

## 🤝 Contributing

We welcome contributions! Please see our [contributing guidelines](CONTRIBUTING.md) for details.

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 👥 Authors

- **Jin Lei** - *Lead Developer*
- **Antonio M. Moro** - *Co-Developer*

## 🔗 Links

- 🌐 **Website**: [smoothie.fewbody.com](https://smoothie.fewbody.com)
<!-- - 📚 **Documentation**: [Full documentation](https://smoothie.fewbody.com/docs) -->
- 🛠️ **Input Generator**: [Web-based input generator](https://smoothie.fewbody.com/generator.html)

---

<div align="center">
Made with ❤️ by the SMOOTHIE team
</div>