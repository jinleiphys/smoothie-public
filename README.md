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

The GUI provides an intuitive interface with forms, real-time plotting, and output monitoring. See the [GUI section](#graphical-user-interface-gui) below for details.

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

## 🖥️ Graphical User Interface (GUI)

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

```bash
# Run the setup script
chmod +x setup_gui.sh
./setup_gui.sh
```

The script will automatically:
1. ✅ Install conda (Miniconda) if not present
2. ✅ Install gfortran compiler if needed
3. ✅ Detect and configure LAPACK/BLAS libraries
4. ✅ Create `smoothie_gui` conda environment
5. ✅ Install Python dependencies (PySide6, matplotlib, numpy)
6. ✅ Build SMOOTHIE Fortran code with optimized settings
7. ✅ Create a launcher script (`run_smoothie_gui.sh`)

### Running the GUI

**Recommended method** - Use the launcher script:
```bash
./run_smoothie_gui.sh
```

**Alternative method** - Manual launch:
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

## 📋 Input File Structure

SMOOTHIE uses Fortran namelist format with the following main sections:

### 🔧 &GLOBAL - Global Parameters

<details>
<summary><strong>Basic Settings</strong></summary>

| Parameter | Description | Default |
|-----------|-------------|---------|
| `hcm` | Radial step size (fm) for wave function integration | 0.05 |
| `rmax` | Maximum radius (fm) for all channels | 50.0 |
| `lmax` | Maximum orbital angular momentum (λₐ, λᵦ) | - |
| `lmin` | Minimum orbital angular momentum | 0 |
| `lxmax` | Maximum lₓ for x+A channel | lmax |
| `elab` | Laboratory energy (MeV) of projectile "a" | - |

</details>

<details>
<summary><strong>Angular Distribution</strong></summary>

| Parameter | Description | Default |
|-----------|-------------|---------|
| `thmin`, `thmax` | Min/max scattering angles (degrees) | - |
| `thinc` | Angular increment (degrees) for output | - |
| `nx` | Gaussian quadrature points for angular integration | 34 |

</details>

<details>
<summary><strong>DWBA Methods</strong></summary>

| Value | Description |
|-------|-------------|
| `dwba=1` | No intrinsic spins, rbx integration variable |
| `dwba=2` | No intrinsic spins, rb integration variable |
| `dwba=3` | With spins, rb integration variable |
| `dwba=4` | Enhanced method with Lagrange mesh |
| `dwba=5` | Lagrange mesh with rb variable |

</details>

### 🎯 &SYSTEM - System Parameters

<details>
<summary><strong>Particle Properties</strong></summary>

| Parameter | Description |
|-----------|-------------|
| `namep`, `namet`, `nameb`, `namex` | Names of projectile, target, detected, and undetected particles |
| `massp`, `masst`, `massb`, `massx` | Masses (amu) of each particle |
| `zp`, `zt`, `zb`, `zx` | Charges of each particle |
| `jp`, `jt`, `jb`, `jx` | Spins of each particle |

</details>

<details>
<summary><strong>Bound State Properties</strong></summary>

| Parameter | Description |
|-----------|-------------|
| `lbx` | Orbital angular momentum of b-x within projectile bound state |
| `sbx` | Total spin coupling (jb + jx) for detected and undetected particles |
| `nodes` | Number of nodes in the b-x bound state wave function |
| `be` | Binding energy (MeV) of the b-x system (positive value) |

</details>

### 📊 &OUTGOING - Output Energy Parameters

| Parameter | Description |
|-----------|-------------|
| `ecmbmin` | Minimum energy (MeV) of outgoing particle b |
| `ecmbmax` | Maximum energy (MeV) of outgoing particle b |
| `ecmbh` | Energy step size (MeV) for output |

### ⚛️ &POTENTIAL - Optical Potentials

<details>
<summary><strong>System Identification (kp1)</strong></summary>

| Value | Description |
|-------|-------------|
| `'a'` | Projectile + target (a+A) channel |
| `'b'` | Detected + recoil (b+B) channel |
| `'t'` | Detected + target (b+A) channel |
| `'x'` | Undetected + target (x+A) channel |
| `'p'` | Detected + undetected (b+x) bound state |

</details>

<details>
<summary><strong>Potential Models (ptype)</strong></summary>

| Value | Description |
|-------|-------------|
| `1` | Woods-Saxon potential with full parameter set |
| `2` | Gaussian potential |
| `3` | Koning-Delaroche (KD02) nucleon-nucleus global potential |
| `4` | CH89 nucleon-nucleus global potential |
| `41-44` | Read potential from external files (fort.41 to fort.44) |

</details>

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

<details>
<summary><strong>View Input File</strong></summary>

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

</details>

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