#!/bin/bash

# SMOOTHIE GUI Launcher Script

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Get conda base path
CONDA_BASE=$(conda info --base 2>/dev/null)
if [ -z "$CONDA_BASE" ]; then
    echo "Error: conda not found!"
    exit 1
fi

# Activate environment
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate smoothie_gui

# Launch GUI
cd "$SCRIPT_DIR/smoothie_gui"
python main.py
