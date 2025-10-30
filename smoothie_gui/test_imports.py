#!/usr/bin/env python3
"""
Quick test to verify all GUI imports work correctly
"""

import sys

def test_imports():
    """Test that all required modules can be imported"""
    
    print("Testing SMOOTHIE GUI imports...\n")
    
    tests = {
        'PySide6': 'Qt framework',
        'matplotlib': 'Plotting library',
        'numpy': 'Numerical computing',
    }
    
    failed = []
    
    for module, description in tests.items():
        try:
            __import__(module)
            print(f"✓ {module:15s} - {description}")
        except ImportError as e:
            print(f"✗ {module:15s} - FAILED: {e}")
            failed.append(module)
    
    print("\nTesting local modules...\n")
    
    local_modules = [
        'main_window',
        'input_panel',
        'plot_widget',
        'log_widget',
        'runner',
        'styles'
    ]
    
    for module in local_modules:
        try:
            __import__(module)
            print(f"✓ {module}")
        except Exception as e:
            print(f"✗ {module} - FAILED: {e}")
            failed.append(module)
    
    print(f"\n{'='*50}")
    if not failed:
        print("✓ All imports successful!")
        print("\nYou can now run: python main.py")
        return 0
    else:
        print(f"✗ {len(failed)} import(s) failed: {', '.join(failed)}")
        print("\nPlease install missing dependencies:")
        print("  pip install -r requirements.txt")
        return 1

if __name__ == '__main__':
    sys.exit(test_imports())
