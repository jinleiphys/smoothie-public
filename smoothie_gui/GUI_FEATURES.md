# SMOOTHIE GUI Features & Design

## Modern Design Philosophy

The SMOOTHIE GUI is built with a modern, user-friendly approach inspired by contemporary scientific software and macOS design principles.

### Visual Design

- **Clean Layout**: Organized split-panel design separating inputs from outputs
- **Tab-Based Organization**: Parameters grouped logically by physics sections
- **Color-Coded Feedback**: Visual status indicators for logs and messages
- **Responsive Interface**: Adapts to different window sizes
- **Professional Typography**: SF Pro Display font family for clarity

### Theme Support

#### Light Theme (Default)
- Clean white backgrounds
- Blue accents (#007AFF)
- High contrast for readability
- Perfect for well-lit environments

#### Dark Theme
- Dark gray backgrounds (#2d2d2d)
- Reduced blue accents (#0A84FF)
- Reduced eye strain
- Ideal for extended work sessions

## Core Components

### 1. Input Panel (Left Side)

**Purpose**: Comprehensive parameter input with validation

**Features**:
- Four organized tabs:
  - **Global**: Mesh, angular momentum, energy, DWBA method
  - **System**: Particle definitions and quantum numbers
  - **Outgoing**: Energy ranges for detected particles
  - **Potentials**: Five potential channels with full parameters

- **Smart Widgets**:
  - Spin boxes with appropriate ranges
  - Combo boxes for categorical selections
  - Real-time validation
  - Tooltips and help text
  - Scroll support for long forms

- **Auto-Population**:
  - Load from existing input files
  - Parse Fortran namelist format
  - Intelligent defaults for each potential type

### 2. Plot Widget (Right Top Tab)

**Purpose**: Real-time visualization of calculation results

**Features**:
- **Multiple Plot Types**:
  1. Cross Section vs Angle (logarithmic scale)
  2. Cross Section vs Energy (linear scale)
  3. Angular Distribution (polar plot)
  4. Energy Spectrum (bar chart)

- **Matplotlib Integration**:
  - Full navigation toolbar (zoom, pan, save)
  - High-quality publication-ready plots
  - Export to PNG, PDF, SVG
  - Customizable appearance

- **Auto-Refresh**:
  - Automatically loads results after calculation
  - Manual refresh button for reloading
  - Handles multiple output file formats

### 3. Log Widget (Right Bottom Tab)

**Purpose**: Real-time monitoring of SMOOTHIE output

**Features**:
- **Color-Coded Messages**:
  - <span style="color: #0066CC">**[INFO]**</span> Blue - Informational messages
  - <span style="color: #FF8800">**[WARNING]**</span> Orange - Warning messages
  - <span style="color: #CC0000">**[ERROR]**</span> Red - Error messages
  - <span style="color: #009900">**[SUCCESS]**</span> Green - Success messages
  - Black - SMOOTHIE stdout/stderr

- **Monospace Font**: Easy to read terminal output
- **Auto-Scroll**: Follows output in real-time
- **Searchable**: Standard text search (Cmd/Ctrl+F)
- **Copy-Paste**: Full text selection and copying

### 4. Menu Bar & Toolbar

**Menu Bar**:
- **File**: New, Open, Save, Save As, Load Example, Exit
- **Run**: Run SMOOTHIE, Stop, Clear Log
- **View**: Theme selection (Light/Dark)
- **Help**: About, Documentation

**Toolbar** (Quick Access):
- New File
- Open File
- Save File
- Load Example
- **Run** (Green button with Ctrl+R)
- **Stop** (Red button with Ctrl+.)

### 5. Status Bar

**Features**:
- Current operation status
- Success/error notifications
- File information
- Persistent across sessions

## Technical Architecture

### Class Hierarchy

```
QMainWindow (main_window.py)
├── InputPanel (input_panel.py)
│   ├── QTabWidget
│   │   ├── GlobalTab (QScrollArea)
│   │   ├── SystemTab (QScrollArea)
│   │   ├── OutgoingTab (QWidget)
│   │   └── PotentialTab (QTabWidget)
│   └── Form Widgets (QSpinBox, QDoubleSpinBox, etc.)
├── PlotWidget (plot_widget.py)
│   ├── QComboBox (plot type selector)
│   ├── FigureCanvas (matplotlib)
│   └── NavigationToolbar (matplotlib)
├── LogWidget (log_widget.py)
│   └── QTextEdit (read-only, styled)
└── SmoothieRunner (runner.py)
    └── QProcess (subprocess management)
```

### Signal-Slot Architecture

```python
# Process signals
SmoothieRunner.output_ready → LogWidget.append_output
SmoothieRunner.error_ready → LogWidget.append_error
SmoothieRunner.started → MainWindow.on_calculation_started
SmoothieRunner.finished → MainWindow.on_calculation_finished
                       → PlotWidget.load_results

# UI signals
MainWindow.run_action.triggered → run_smoothie()
MainWindow.save_action.triggered → save_file()
PlotWidget.plot_type.changed → update_plot()
```

### File I/O

**Input File Parsing**:
```python
# Regex-based namelist parsing
&GLOBAL\s+(.*?)\s*/  # Captures GLOBAL parameters
&SYSTEM\s+(.*?)\s*/  # Captures SYSTEM parameters
# etc.
```

**Output File Generation**:
```python
# Fortran namelist format
"""
NAMELIST
&GLOBAL hcm=0.05 lmax=25 ... /
&SYSTEM namep='d' massp=2.0 ... /
...
"""
```

## User Workflow

### Typical Session

1. **Launch** GUI → `./run_smoothie_gui.sh`
2. **Load** example or create new input
3. **Review/Modify** parameters in tabs
4. **Save** input file (optional)
5. **Run** calculation (Ctrl+R)
6. **Monitor** progress in Output Log
7. **Analyze** results in Plot tab
8. **Export** plots for publication

### Advanced Workflow

1. **Parameter Study**:
   - Load base input
   - Modify energy/angle ranges
   - Run multiple calculations
   - Compare results

2. **Potential Optimization**:
   - Try different potential models
   - Adjust depths/radii
   - Compare with data
   - Iterate to best fit

3. **Batch Processing**:
   - Save multiple input files
   - Run sequentially
   - Collect output files
   - Post-process results

## Performance Optimizations

### GUI Responsiveness

- **Threaded Execution**: SMOOTHIE runs in separate QProcess
- **Non-Blocking UI**: GUI remains responsive during calculation
- **Chunked Output**: Log updates in real-time without freezing
- **Lazy Loading**: Plots only refresh when tab is visible

### Memory Management

- **Efficient Widgets**: Reuse existing widgets vs. recreation
- **Plot Caching**: Results cached until new calculation
- **Bounded Logs**: Optional log size limits for long runs

## Accessibility Features

- **Keyboard Navigation**: Full keyboard support (Tab, Arrow keys)
- **Shortcuts**: Intuitive Ctrl+Key combinations
- **Clear Labels**: All inputs properly labeled
- **Help Text**: Tooltips and descriptions
- **High Contrast**: Both themes meet WCAG standards
- **Readable Fonts**: Large, clear typography

## Cross-Platform Support

### macOS
- Native look and feel with SF Pro Display
- Qt6 backend for modern macOS
- Retina display support
- Cmd key shortcuts

### Linux
- GTK/Qt theme integration
- Font fallbacks (Segoe UI, etc.)
- Ctrl key shortcuts
- Wayland/X11 compatible

### Windows
- Windows 10/11 modern styling
- DirectWrite font rendering
- Windows shortcuts
- High DPI scaling

## Future Enhancements

Potential features for future versions:

- [ ] **Batch Mode**: Queue multiple calculations
- [ ] **Data Import**: Load experimental data for comparison
- [ ] **Fit Module**: Automatic parameter optimization
- [ ] **3D Plots**: Surface plots for parameter scans
- [ ] **Export Presets**: Save favorite parameter sets
- [ ] **Undo/Redo**: Parameter change history
- [ ] **Diff Viewer**: Compare input files
- [ ] **Remote Execution**: Run on HPC clusters
- [ ] **Plugin System**: Custom potential models
- [ ] **Jupyter Integration**: Export to notebooks

## Dependencies

### Core Requirements
```
PySide6 >= 6.5.0      # Qt6 Python bindings
matplotlib >= 3.7.0   # Scientific plotting
numpy >= 1.24.0       # Numerical arrays
```

### Optional Enhancements
```
scipy >= 1.10.0       # Advanced numerical operations
pandas >= 2.0.0       # Data analysis (future)
```

### Build Requirements
```
gfortran >= 4.8       # Fortran compiler
LAPACK/BLAS          # Linear algebra libraries
make                 # Build system
```

## Code Quality

### Standards
- **PEP 8**: Python style guide compliance
- **Type Hints**: Modern Python typing (3.8+)
- **Docstrings**: All classes and methods documented
- **Comments**: Complex logic explained inline

### Testing Strategy
- Manual testing on macOS, Linux, Windows
- Example calculations verified against reference
- UI responsiveness tested on various screen sizes
- Theme consistency checked in both modes

## License & Credits

**License**: MIT (same as SMOOTHIE)

**Built With**:
- Qt for Python (PySide6)
- Matplotlib for plotting
- NumPy for numerical operations

**Inspired By**:
- Modern scientific software UIs
- Apple Human Interface Guidelines
- Material Design principles
- Web-based input generator

---

*For questions, issues, or contributions, please visit the SMOOTHIE repository or website.*
