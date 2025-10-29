# Working Directory Feature

## Overview

The GUI now asks you to select a working directory before running each SMOOTHIE calculation. This allows you to:

- Organize different calculations in separate folders
- Keep input files and output files together
- Easily compare results from different runs
- Avoid cluttering a single directory with multiple outputs

## How It Works

### When You Click "Run"

1. **Directory Selection Dialog** appears
   - Choose an existing folder OR create a new one
   - Default location: `~/Desktop/code/smoothie-public/smoothie/test`
   - The GUI remembers your last selection for convenience

2. **Files Created in That Directory**:
   - `smoothie_input.in` - Your input parameters
   - `fort.20` - Cross section results
   - `fort.21` - Angular distributions
   - `fort.22` - Energy spectra
   - `fort.1-4` - Standard SMOOTHIE outputs
   - `fort.910-919` - Detailed channel data

3. **Automatic Result Loading**:
   - After calculation completes, GUI automatically looks in that directory
   - Plots are loaded from the output files
   - "Refresh" button reloads from the same directory

## Usage Example

### Organize Calculations by Project

```
~/calculations/
  ├── deuteron_25MeV/
  │   ├── smoothie_input.in
  │   ├── fort.20
  │   ├── fort.21
  │   └── ...
  │
  ├── deuteron_50MeV/
  │   ├── smoothie_input.in
  │   ├── fort.20
  │   └── ...
  │
  └── beryllium_30MeV/
      ├── smoothie_input.in
      └── ...
```

### Workflow

1. **First Calculation**:
   ```
   - Click Run
   - Select/Create: ~/calculations/deuteron_25MeV/
   - Calculation runs, files saved there
   - Results plotted
   ```

2. **Second Calculation**:
   ```
   - Modify parameters in GUI
   - Click Run
   - Select/Create: ~/calculations/deuteron_50MeV/
   - New calculation runs in separate folder
   - Both results preserved
   ```

3. **Compare Results**:
   - Each folder has complete input + output
   - Easy to replot or reanalyze later
   - No confusion about which parameters produced which results

## Benefits

### Organization
- Clean separation of different calculations
- Easy to find results later
- Can add notes/plots in each folder

### Reproducibility
- Input file saved with results
- Complete record of calculation
- Easy to rerun with same parameters

### Flexibility
- Run multiple parameter studies
- Compare different reactions
- Archive old calculations

## Tips

### Create Descriptive Folder Names

Good examples:
- `d_Nb93_25MeV_DWBA1`
- `Be11_Zn64_test1`
- `parameter_scan_energy`

### Use Subfolders for Projects

```
~/smoothie_calculations/
  ├── project1_deuteron_breakup/
  │   ├── energy_scan/
  │   │   ├── 20MeV/
  │   │   ├── 25MeV/
  │   │   └── 30MeV/
  │   └── angle_scan/
  │       ├── fine_mesh/
  │       └── coarse_mesh/
  └── project2_beryllium/
```

### Quick Access
- The GUI remembers your last directory
- Press Ctrl+R to run again in same location
- Create shortcuts to frequently used directories

## Technical Details

### What Changed

1. **main_window.py**:
   - Added `self.working_directory` to store current directory
   - Modified `run_smoothie()` to show directory selection dialog
   - Input file saved as `smoothie_input.in` in chosen directory
   - Working directory passed to runner

2. **runner.py**:
   - Updated `run()` to accept `working_dir` parameter
   - Process runs in specified directory
   - Output files created there automatically

3. **plot_widget.py**:
   - Added `self.working_directory` to remember location
   - Updated `load_results()` to accept directory parameter
   - Looks for fort.* files in specified directory
   - Refresh button uses remembered directory

### Backward Compatibility

The changes are backward compatible:
- If no directory specified, uses current directory
- Old behavior preserved for programmatic use
- All parameters optional with sensible defaults

## Cancelling

If you click Run but then cancel the directory selection:
- Calculation does not start
- No error message (normal behavior)
- Click Run again when ready

## Future Enhancements

Potential improvements:
- [ ] Project management system
- [ ] Recently used directories list
- [ ] Batch mode: queue multiple directories
- [ ] Compare results from different folders
- [ ] Export summary across multiple calculations

---

**Added in**: Version 1.1
**Last Updated**: 2025-10-30
