# SMOOTHIE GUI - Successfully Installed! 🎉

The SMOOTHIE GUI has been successfully installed and tested on your system.

## ✅ Installation Complete

All components are working:
- ✅ PySide6 (Qt framework) installed
- ✅ matplotlib (plotting) installed  
- ✅ numpy (numerical computing) installed
- ✅ GUI modules loaded successfully
- ✅ Application launches without errors

## 🚀 Launch the GUI

Simply run:

```bash
cd /Users/jinlei/Desktop/code/smoothie-public
./run_smoothie_gui.sh
```

Or manually:

```bash
conda activate smoothie_gui
cd smoothie_gui
python main.py
```

## 📝 Quick Usage

1. **Load Example**:
   - Click `File → Load Example` (or toolbar button)
   - This loads the test.in deuteron breakup example

2. **Review/Modify Parameters**:
   - Navigate through tabs: Global, System, Outgoing, Potentials
   - Adjust values as needed

3. **Run Calculation**:
   - Click `Run` button (or press Ctrl+R)
   - Watch output in "Output Log" tab
   - Results auto-display in "Plot" tab

4. **Save Your Work**:
   - Click `File → Save` (Ctrl+S) to save input file
   - Export plots using matplotlib toolbar

## 🎨 Features Available

### Input Panel (Left)
- **Global**: Energy, DWBA method, angular/radial mesh
- **System**: Projectile, target, particles, quantum numbers
- **Outgoing**: Energy ranges
- **Potentials**: 5 channels with full parameter control

### Output Panel (Right)
- **Plot Tab**: 4 plot types (cross section, angular dist, etc.)
- **Log Tab**: Real-time color-coded output

### Themes
- Light theme (default): `View → Theme → Light`
- Dark theme: `View → Theme → Dark`

## ⚠️ Font Note

You may see this Qt warning when launching:
```
qt.qpa.fonts: Populating font family aliases took XXX ms. 
Replace uses of missing font family "SF Pro Display"...
```

**This is harmless!** Qt automatically uses your system's default font. The GUI works perfectly. If you want to eliminate the warning, you can install the SF Pro Display font from Apple, or just ignore it.

## 📚 Documentation

- **Quick Start**: `GUI_QUICKSTART.md`
- **Full Guide**: `README.md`
- **Features**: `GUI_FEATURES.md`
- **Visual Guide**: `VISUAL_GUIDE.md`

## 🎯 Example Workflow

1. Launch GUI
2. Load example (File menu or toolbar)
3. Click Run (or Ctrl+R)
4. Switch to Plot tab when done
5. Choose plot type from dropdown
6. Save/export as needed

## 💡 Tips

- **Keyboard Shortcuts**:
  - `Ctrl+R`: Run calculation
  - `Ctrl+S`: Save input
  - `Ctrl+O`: Open file
  - `Ctrl+N`: New file
  - `Ctrl+.`: Stop calculation

- **Plot Tools**:
  - Use matplotlib toolbar to zoom/pan
  - Save plots as PNG, PDF, or SVG
  - Multiple plot types available

- **Themes**:
  - Try dark theme for extended work sessions
  - Switch via View → Theme menu

## 🔧 Troubleshooting

If you encounter any issues:

1. **GUI won't start**: Check Python version is 3.8+
   ```bash
   conda activate smoothie_gui
   python --version
   ```

2. **SMOOTHIE not found**: Build it first
   ```bash
   cd smoothie
   make clean && make
   ```

3. **Import errors**: Reinstall dependencies
   ```bash
   pip install --upgrade -r requirements.txt
   ```

4. **Font warnings**: Ignore them or install SF Pro Display font

## 🎉 Enjoy!

You now have a beautiful, modern GUI for SMOOTHIE calculations. Explore the features, try different reactions, and enjoy the integrated plotting and monitoring capabilities!

For questions or issues, see the documentation files or visit:
https://smoothie.fewbody.com

---

**Built with**: PySide6 • matplotlib • Python  
**License**: MIT (same as SMOOTHIE)
