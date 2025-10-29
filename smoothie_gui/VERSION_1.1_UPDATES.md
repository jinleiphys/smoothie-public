# SMOOTHIE GUI Version 1.1 - Update Summary

## What's New

### 🎯 Working Directory Selection Feature

**Before**: Calculations ran in a temporary location, files scattered

**Now**: When you click "Run", you select where to save everything!

#### How It Works:

1. Click **Run** button (or press Ctrl+R)
2. **Dialog appears**: "Select Working Directory for Calculation"
3. Choose or create a folder (e.g., `~/calculations/deuteron_25MeV`)
4. Everything saved there:
   - `smoothie_input.in` - Your input parameters
   - `fort.20, fort.21, fort.22` - Results
   - All other SMOOTHIE output files

5. Results automatically loaded from that directory
6. Next time you run, starts in same location (convenient!)

#### Benefits:

✅ **Organized**: Each calculation in its own folder
✅ **Reproducible**: Input + output together
✅ **Comparable**: Easy to compare different runs
✅ **Clean**: No file clutter
✅ **Professional**: Archive and share complete calculations

### 📝 Updated Window Title

**New title**: `SMOOTHIE - IAV Nonelastic breakup calculations`

More descriptive and accurate for the physics being calculated (Ichimura-Austern-Vincent nonelastic breakup reactions).

### ℹ️ Enhanced About Dialog

Updated "About" dialog (Help → About) with:
- IAV Nonelastic Breakup Calculations subtitle
- Version 1.1 notation
- Improved formatting

## Example Workflow

### Scenario: Parameter Study

```bash
# Run 1: 25 MeV energy
- Load example or set parameters
- Click Run
- Select: ~/calculations/d_Nb_25MeV/
- Results saved there

# Run 2: 30 MeV energy
- Change energy to 30 MeV
- Click Run
- Select: ~/calculations/d_Nb_30MeV/
- Results saved there (Run 1 preserved!)

# Run 3: 35 MeV energy
- Change energy to 35 MeV
- Click Run
- Select: ~/calculations/d_Nb_35MeV/
- Results saved there

# Now you have:
~/calculations/
  ├── d_Nb_25MeV/
  │   ├── smoothie_input.in
  │   ├── fort.20
  │   └── ...
  ├── d_Nb_30MeV/
  │   ├── smoothie_input.in
  │   ├── fort.20
  │   └── ...
  └── d_Nb_35MeV/
      ├── smoothie_input.in
      ├── fort.20
      └── ...
```

Perfect for parameter scans, comparison studies, or publications!

## Technical Changes

### Modified Files:

1. **main_window.py**:
   - Added `self.working_directory` instance variable
   - Updated `run_smoothie()` to show directory selection dialog
   - Changed window title to include "IAV Nonelastic breakup"
   - Updated About dialog with version 1.1
   - Input file now saved as `smoothie_input.in` in chosen directory

2. **runner.py**:
   - Modified `run()` method to accept `working_dir` parameter
   - Process now runs in specified directory
   - Backward compatible (optional parameter)

3. **plot_widget.py**:
   - Added `self.working_directory` instance variable
   - Updated `load_results()` to accept working directory
   - Looks for output files in specified directory
   - Refresh button uses remembered directory

### Backward Compatibility:

✅ All changes are backward compatible
✅ Optional parameters with sensible defaults
✅ Previous behavior preserved when directory not specified

## Usage Tips

### Organize by Project

```
~/smoothie_calculations/
  ├── thesis_chapter3/
  │   ├── deuteron_tests/
  │   └── beryllium_tests/
  ├── paper_2025/
  │   ├── figure1_data/
  │   └── figure2_data/
  └── preliminary/
```

### Descriptive Names

Good folder names tell the story:
- `d_Pb208_56MeV_DWBA1`
- `Be11_Zn64_angular_scan`
- `test_new_potential_v2`

### Quick Rerun

The GUI remembers your last directory:
- Modify parameters
- Click Run (Ctrl+R)
- Same folder selected by default
- Easy to overwrite or choose new location

## Migration from v1.0

No migration needed! Just start using the new feature:

1. **First run after update**:
   - Click Run → Select directory → Works!

2. **Existing calculations**:
   - Still accessible in their original locations
   - No files moved or changed

3. **Learning curve**:
   - Zero! Just choose a folder when prompted
   - Very intuitive interface

## Troubleshooting

### Q: Can I cancel the directory selection?

**A**: Yes! Click Cancel and nothing happens. Click Run again when ready.

### Q: Can I use the same directory twice?

**A**: Yes! Files will be overwritten. Consider adding `_v2` or `_test2` to folder name.

### Q: Where should I save calculations?

**A**: Anywhere you want! Suggestions:
- `~/calculations/` - Simple, organized
- `~/Documents/SMOOTHIE/` - With other docs
- `~/Desktop/current_project/` - Easy access
- External drive for large studies

### Q: Do I need to select directory every time?

**A**: Yes, but it remembers your last choice as the default, so it's quick!

### Q: Can I still use the old behavior?

**A**: The GUI always asks for a directory now. This is better for organization!

## What's Next?

Potential future features:
- Recent directories list (quick access)
- Project management (collections of calculations)
- Batch mode (queue multiple runs)
- Compare results across directories
- Export summary reports

## Version History

### Version 1.1 (2025-10-30)
- ✨ NEW: Working directory selection
- 📝 Updated window title to "IAV Nonelastic breakup calculations"
- ℹ️ Enhanced About dialog
- 🐛 Fixed initialization order bug

### Version 1.0 (2025-10-29)
- 🎉 Initial release
- Complete GUI with input forms, plotting, and logging
- Light and dark themes
- File operations (load, save, examples)
- Real-time output monitoring

---

**Enjoy the improved workflow!** 🚀

For questions or feedback, see README.md or visit https://smoothie.fewbody.com
