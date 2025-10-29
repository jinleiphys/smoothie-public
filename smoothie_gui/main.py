#!/usr/bin/env python3
"""
SMOOTHIE GUI - Modern interface for SMOOTHIE quantum scattering calculations
"""

import sys
from PySide6.QtWidgets import QApplication
from PySide6.QtCore import Qt
from main_window import MainWindow


def main():
    """Main entry point for the SMOOTHIE GUI application"""
    # Enable high DPI scaling
    QApplication.setHighDpiScaleFactorRoundingPolicy(
        Qt.HighDpiScaleFactorRoundingPolicy.PassThrough
    )

    app = QApplication(sys.argv)
    app.setApplicationName("SMOOTHIE")
    app.setOrganizationName("SMOOTHIE Team")

    # Set application-wide font
    from PySide6.QtGui import QFont
    font = QFont("SF Pro Display", 10)
    app.setFont(font)

    window = MainWindow()
    window.show()

    sys.exit(app.exec())


if __name__ == "__main__":
    main()
