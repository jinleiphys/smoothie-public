# SMOOTHIE GUI

Modern graphical user interface for SMOOTHIE quantum scattering calculations.

## Features

- **Input Panel**: Comprehensive form-based input with tabs for different parameter sections (Global, System, Outgoing, Potentials)
- **Plotting Area**: Real-time visualization of calculation results with multiple plot types
- **Output Log**: Color-coded log display with real-time monitoring of SMOOTHIE output
- **Modern Design**: Beautiful, responsive interface with light and dark theme support
- **File Operations**: Load/save input files, load examples, export results
- **Integrated Execution**: Run SMOOTHIE directly from the GUI with real-time feedback

## Screenshots

The GUI features a split-panel layout:
- Left panel: Parameter input forms organized in tabs
- Right panel: Results plotting and output logs

## Installation

### Quick Setup (Recommended)

Run the automated setup script:

```bash
cd /path/to/smoothie-public
chmod +x setup_gui.sh
./setup_gui.sh
```

This script will:
1. Create a conda environment named `smoothie_gui`
2. Install all Python dependencies
3. Build the SMOOTHIE Fortran code
4. Set up the environment for immediate use

### Manual Setup

1. **Create conda environment**:
```bash
conda create -n smoothie_gui python=3.10
conda activate smoothie_gui
```

2. **Install Python dependencies**:
```bash
cd smoothie_gui
pip install -r requirements.txt
```

3. **Build SMOOTHIE**:
```bash
cd ../smoothie
cp ../make.inc.gfortran ../make.inc  # or use make.inc.ifort for Intel compiler
make clean
make
```

## Usage

### Running the GUI

```bash
# Activate the environment
conda activate smoothie_gui

# Navigate to GUI directory
cd smoothie_gui

# Launch the application
python main.py
```

### Using the Application

1. **Input Parameters**:
   - Navigate through the tabs (Global, System, Outgoing, Potentials)
   - Fill in the required parameters for your calculation
   - Or load an example using File → Load Example

2. **Run Calculation**:
   - Click the "Run" button in the toolbar or press Ctrl+R
   - Monitor progress in the Output Log tab
   - Results will automatically appear in the Plot tab when complete

3. **View Results**:
   - Switch to the Plot tab to visualize results
   - Choose different plot types from the dropdown menu
   - Use the matplotlib toolbar to zoom, pan, and save plots

4. **Save/Load**:
   - Save current inputs: File → Save (Ctrl+S)
   - Load previous inputs: File → Open (Ctrl+O)
   - Export plots using the matplotlib toolbar

### Keyboard Shortcuts

- `Ctrl+N`: New file
- `Ctrl+O`: Open file
- `Ctrl+S`: Save file
- `Ctrl+R`: Run SMOOTHIE
- `Ctrl+.`: Stop calculation
- `Ctrl+Q`: Quit application

## Themes

The GUI supports both light and dark themes:
- Switch themes via View → Theme menu
- Light theme: Clean, bright interface (default)
- Dark theme: Easy on the eyes for extended use

## Requirements

- Python 3.8 or higher
- PySide6 (Qt for Python)
- matplotlib
- numpy
- SMOOTHIE compiled executable

## Troubleshooting

### SMOOTHIE executable not found

If you see an error about the SMOOTHIE executable:
1. Ensure SMOOTHIE is built: `cd smoothie && make`
2. Check the executable path in `main_window.py` (line ~275)
3. Update the path if your installation differs

### Import errors

If you encounter import errors:
```bash
conda activate smoothie_gui
pip install --upgrade -r requirements.txt
```

### GUI doesn't start

On macOS, you may need to install system dependencies:
```bash
brew install qt6
```

## Development

### Project Structure

```
smoothie_gui/
├── main.py              # Application entry point
├── main_window.py       # Main window and menu/toolbar
├── input_panel.py       # Parameter input forms
├── plot_widget.py       # Matplotlib plotting widget
├── log_widget.py        # Output log display
├── runner.py            # SMOOTHIE execution handler
├── styles.py            # Theme definitions
├── requirements.txt     # Python dependencies
└── README.md           # This file
```

### Contributing

To contribute to the GUI:
1. Follow PEP 8 style guidelines
2. Use Qt's signal/slot mechanism for inter-widget communication
3. Maintain compatibility with both light and dark themes
4. Test on multiple platforms when possible

## License

This GUI is part of the SMOOTHIE project and follows the same MIT license.

## Support

For issues or questions:
- Check the main SMOOTHIE documentation
- Visit https://smoothie.fewbody.com
- Open an issue on GitHub

## Credits

Built with:
- PySide6 (Qt for Python)
- matplotlib for scientific plotting
- Modern design inspired by macOS and web interfaces
