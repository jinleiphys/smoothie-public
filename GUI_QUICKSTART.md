# SMOOTHIE GUI Quick Start Guide

This guide will help you get the SMOOTHIE GUI up and running in minutes.

## One-Command Setup

From the `smoothie-public` directory, run:

```bash
./setup_gui.sh
```

This single script will:
1. Create a conda environment with Python 3.10
2. Install all required Python packages (PySide6, matplotlib, numpy)
3. Build the SMOOTHIE Fortran code and all modules
4. Create a launcher script for easy access

## Running the GUI

After setup completes, launch the GUI with:

```bash
./run_smoothie_gui.sh
```

Or manually:

```bash
conda activate smoothie_gui
cd smoothie_gui
python main.py
```

## First Calculation

1. **Load an Example**:
   - Click `File → Load Example` (or the "Load Example" button in toolbar)
   - This loads the default deuteron breakup example

2. **Review Parameters**:
   - Navigate through the tabs: Global, System, Outgoing, Potentials
   - See the preset values for a d+93Nb reaction at 25.5 MeV

3. **Run Calculation**:
   - Click the "Run" button in the toolbar (or press `Ctrl+R`)
   - Watch the output in the "Output Log" tab
   - Progress updates appear in real-time

4. **View Results**:
   - Switch to the "Plot" tab
   - Results automatically display when calculation completes
   - Use dropdown to switch between different plot types:
     - Cross Section vs Angle
     - Cross Section vs Energy
     - Angular Distribution
     - Energy Spectrum

## GUI Overview

### Layout

```
┌─────────────────────────────────────────────────────────────┐
│ File  Run  View  Help                    [Toolbar Buttons]  │
├────────────────┬────────────────────────────────────────────┤
│                │  ┌──────────────────────────────────────┐  │
│  Input Panel   │  │ Plot / Output Log (Tabs)             │  │
│  ┌──────────┐  │  │                                      │  │
│  │ Global   │  │  │  [Matplotlib Canvas]                │  │
│  │ System   │  │  │                                      │  │
│  │ Outgoing │  │  │  or                                  │  │
│  │ Potentials│ │  │                                      │  │
│  └──────────┘  │  │  [Color-coded Log Output]            │  │
│                │  └──────────────────────────────────────┘  │
├────────────────┴────────────────────────────────────────────┤
│ Status: Ready                                                │
└─────────────────────────────────────────────────────────────┘
```

### Input Panel (Left Side)

Four tabs organize all input parameters:

- **Global**: Radial steps, angular momentum ranges, energy, DWBA method, angular range, quadrature points
- **System**: Projectile, target, detected/undetected particles, binding energy, quantum numbers
- **Outgoing**: Energy range for outgoing particles
- **Potentials**: Five potential definitions (projectile-target, detected-recoil, etc.)

### Output/Plot Panel (Right Side)

Two tabs for results:

- **Plot**: Real-time visualization with matplotlib
  - Multiple plot types via dropdown
  - Full matplotlib toolbar (zoom, pan, save)
  - Auto-refresh when calculation completes

- **Output Log**: Color-coded terminal output
  - Blue: Info messages
  - Orange: Warnings
  - Red: Errors
  - Green: Success
  - Black: SMOOTHIE output

## Keyboard Shortcuts

| Shortcut | Action |
|----------|--------|
| `Ctrl+N` | New file (reset inputs) |
| `Ctrl+O` | Open input file |
| `Ctrl+S` | Save input file |
| `Ctrl+Shift+S` | Save as... |
| `Ctrl+R` | Run SMOOTHIE |
| `Ctrl+.` | Stop calculation |
| `Ctrl+Q` | Quit application |

## Themes

Switch between light and dark themes via `View → Theme`:

- **Light**: Clean, bright interface (default)
- **Dark**: Reduced eye strain for extended sessions

## File Operations

### Saving Inputs

1. `File → Save` or `Ctrl+S`
2. Choose location and filename (`.in` extension recommended)
3. File is saved in SMOOTHIE namelist format

### Loading Inputs

1. `File → Open` or `Ctrl+O`
2. Select an `.in` file
3. All parameters are populated automatically

### Input File Format

Files are saved in standard SMOOTHIE Fortran namelist format:

```fortran
NAMELIST
&GLOBAL hcm=0.05 lmax=25 elab=25.5 ... /
&SYSTEM namep='d' massp=2.0 ... /
&OUTGOING ecmbmin=2.0 ... /
&POTENTIAL kp1='a' ptype=1 ... /
...
```

## Example Calculations

The GUI includes several example scenarios:

### Deuteron Breakup (Default Example)

Reaction: d + 93Nb → p + n + 93Nb at 25.5 MeV

- DWBA method: 1 (no spins, rbx variable)
- Potentials: Woods-Saxon for d+Nb, CH89 global for p+Nb and n+Nb
- Binding form factor: Gaussian

**To run**: Load Example → Run

### Custom Calculations

1. Start with an example as template
2. Modify parameters as needed:
   - Change projectile/target
   - Adjust energy
   - Select different potentials
   - Modify angular/energy ranges
3. Save your custom input
4. Run and analyze

## Troubleshooting

### "SMOOTHIE executable not found"

**Solution**: Build SMOOTHIE first:
```bash
cd smoothie
make clean
make
```

### GUI won't start

**Check Python version**:
```bash
conda activate smoothie_gui
python --version  # Should be 3.8+
```

**Reinstall dependencies**:
```bash
pip install --upgrade -r smoothie_gui/requirements.txt
```

### Calculation errors

1. Check Output Log tab for specific error messages
2. Verify input parameters are physical
3. Try loading the default example to ensure SMOOTHIE works
4. Check that all required files exist in working directory

### Display issues on macOS

If Qt rendering issues occur:
```bash
brew install qt6
```

## Advanced Features

### DWBA Methods

Choose calculation method in Global tab:

1. **DWBA=1**: No spins, rbx variable (fastest, most common, PLM optimized)
2. **DWBA=2**: No spins, rb variable
3. **DWBA=3**: With spins, rb variable (includes spin degrees of freedom)
4. **DWBA=4**: Enhanced Lagrange mesh method
5. **DWBA=5**: Lagrange mesh with rb variable

### Potential Models

Each of 5 potential channels can use:

- **Woods-Saxon**: Traditional optical potential
- **Gaussian**: Smooth form factor for bound states
- **KD02 Global**: Koning-Delaroche global optical potential
- **CH89 Global**: Chapel Hill 89 global potential
- **Global Deuteron**: Specialized deuteron potential
- **Read from fort.XX**: Custom tabulated potentials

### Output Files

SMOOTHIE generates multiple output files:

- `fort.20`: Cross sections vs angle
- `fort.21`: Angular distributions
- `fort.22`: Energy spectra
- `fort.910-919`: Detailed channel information

The GUI automatically reads and plots these when available.

## Getting Help

- **GUI Documentation**: `smoothie_gui/README.md`
- **SMOOTHIE Documentation**: `CLAUDE.md` in project root
- **Website**: https://smoothie.fewbody.com
- **Issues**: GitHub repository

## Next Steps

1. **Run the default example** to verify everything works
2. **Explore the interface** and different plot types
3. **Modify parameters** to understand their effects
4. **Try different reactions** by changing projectile/target
5. **Experiment with DWBA methods** and potential models
6. **Compare with experimental data** (load your own data)

Enjoy using SMOOTHIE GUI for your quantum scattering calculations!
