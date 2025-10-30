# Intelligent Path Detection

The SMOOTHIE GUI now includes intelligent path detection to automatically find executables and directories, regardless of the user's installation location.

## Features

### 1. Automatic Repository Detection
The GUI automatically detects the SMOOTHIE repository root by searching upward from the GUI script location. This works across different users and installation paths.

### 2. Smart Executable Discovery
The GUI searches for executables (smoothie and cm2lab) in the following order:

1. **Environment variables**: `SMOOTHIE_EXE` or `CM2LAB_EXE`
2. **Relative to detected repo root**: `<repo>/smoothie/smoothie` or `<repo>/cm2lab/cm2lab`
3. **System PATH**: Any location in your PATH environment variable
4. **Common installation locations**: `/usr/local/bin`, `/opt/smoothie/bin`, `~/bin`

### 3. Dynamic Default Directories
File dialogs automatically use the detected test directory (`<repo>/smoothie/test`) as the default location, making it easier to work with example files.

## Configuration

### Using Environment Variables
If your executables are in a non-standard location, you can set environment variables:

```bash
# For SMOOTHIE
export SMOOTHIE_EXE=/path/to/your/smoothie

# For CM2LAB
export CM2LAB_EXE=/path/to/your/cm2lab

# Then run the GUI
./run_smoothie_gui.sh
```

### For Development
If you're running the GUI from a different location or have multiple SMOOTHIE installations, the GUI will:
- Automatically find the nearest repository
- Display the detected executable paths in the log
- Provide helpful build instructions if executables aren't found

## Benefits

- **No hardcoded paths**: Works for any user without modification
- **Cross-platform**: Adapts to different directory structures
- **Better error messages**: Provides specific instructions when executables aren't found
- **Flexible configuration**: Supports custom installations via environment variables

## Error Messages

If the GUI can't find an executable, it will show:
- The expected location based on the detected repository
- Instructions to build the executable
- Alternative configuration via environment variables
