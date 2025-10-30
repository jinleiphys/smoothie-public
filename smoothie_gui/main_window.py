"""
Main window for SMOOTHIE GUI application
"""

from PySide6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QSplitter,
    QTabWidget, QToolBar, QStatusBar, QFileDialog, QMessageBox, QToolButton
)
from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QAction, QIcon, QKeySequence
import os

from input_panel import InputPanel
from plot_widget import PlotWidget
from log_widget import LogWidget
from runner import SmoothieRunner
from styles import apply_modern_style


class MainWindow(QMainWindow):
    """Main application window with modern layout"""

    def __init__(self):
        super().__init__()
        self.smoothie_runner = SmoothieRunner()
        self.current_file = None
        self.working_directory = None  # Store current working directory
        self.is_running_smoothie = False  # Track whether SMOOTHIE or CM2LAB is running

        self.init_ui()
        self.setup_connections()

    def init_ui(self):
        """Initialize the user interface"""
        self.setWindowTitle("SMOOTHIE - IAV nonelastic breakup calculations")
        self.setGeometry(100, 100, 1600, 900)

        # Apply modern styling (default: light theme, change to "dark" for dark theme)
        apply_modern_style(self)

        # Create central widget with splitter layout FIRST
        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        main_layout = QHBoxLayout(central_widget)
        main_layout.setContentsMargins(0, 0, 0, 0)

        # Create main splitter (horizontal)
        splitter = QSplitter(Qt.Horizontal)

        # Left panel: Input forms
        self.input_panel = InputPanel()
        splitter.addWidget(self.input_panel)

        # Right panel: Tabs for plotting and output
        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)
        right_layout.setContentsMargins(0, 0, 0, 0)

        self.right_tabs = QTabWidget()
        self.right_tabs.setDocumentMode(True)

        # Plot tab
        self.plot_widget = PlotWidget()
        self.right_tabs.addTab(self.plot_widget, "Plot")

        # Log tab
        self.log_widget = LogWidget()
        self.right_tabs.addTab(self.log_widget, "Output Log")

        right_layout.addWidget(self.right_tabs)
        splitter.addWidget(right_widget)

        # Set initial splitter sizes (40% input, 60% output/plot)
        splitter.setSizes([640, 960])

        main_layout.addWidget(splitter)

        # Create status bar
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready")

        # NOW create menu bar and toolbar (after widgets exist)
        self.create_menus()
        self.create_toolbar()

    def create_menus(self):
        """Create menu bar with all menu items"""
        menubar = self.menuBar()

        # File menu
        file_menu = menubar.addMenu("&File")

        new_action = QAction("&New", self)
        new_action.setShortcut(QKeySequence.New)
        new_action.triggered.connect(self.new_file)
        file_menu.addAction(new_action)

        open_action = QAction("&Open...", self)
        open_action.setShortcut(QKeySequence.Open)
        open_action.triggered.connect(self.open_file)
        file_menu.addAction(open_action)

        self.save_action = QAction("&Save", self)
        self.save_action.setShortcut(QKeySequence.Save)
        self.save_action.triggered.connect(self.save_file)
        file_menu.addAction(self.save_action)

        save_as_action = QAction("Save &As...", self)
        save_as_action.setShortcut(QKeySequence.SaveAs)
        save_as_action.triggered.connect(self.save_file_as)
        file_menu.addAction(save_as_action)

        file_menu.addSeparator()

        exit_action = QAction("E&xit", self)
        exit_action.setShortcut(QKeySequence.Quit)
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

        # Run menu
        run_menu = menubar.addMenu("&Run")

        self.run_action = QAction("&Run SMOOTHIE", self)
        self.run_action.setShortcut("Ctrl+R")
        self.run_action.triggered.connect(self.run_smoothie)
        run_menu.addAction(self.run_action)

        self.stop_action = QAction("&Stop", self)
        self.stop_action.setShortcut("Ctrl+.")
        self.stop_action.setEnabled(False)
        self.stop_action.triggered.connect(self.stop_smoothie)
        run_menu.addAction(self.stop_action)

        run_menu.addSeparator()

        clear_log_action = QAction("&Clear Log", self)
        clear_log_action.triggered.connect(self.log_widget.clear)
        run_menu.addAction(clear_log_action)

        # View menu
        view_menu = menubar.addMenu("&View")

        theme_menu = view_menu.addMenu("&Theme")

        light_theme_action = QAction("&Light", self)
        light_theme_action.triggered.connect(lambda: self.change_theme("light"))
        theme_menu.addAction(light_theme_action)

        dark_theme_action = QAction("&Dark", self)
        dark_theme_action.triggered.connect(lambda: self.change_theme("dark"))
        theme_menu.addAction(dark_theme_action)

        # Help menu
        help_menu = menubar.addMenu("&Help")

        about_action = QAction("&About", self)
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)

        docs_action = QAction("&Documentation", self)
        docs_action.triggered.connect(self.show_docs)
        help_menu.addAction(docs_action)

    def create_toolbar(self):
        """Create toolbar with quick action buttons"""
        toolbar = QToolBar("Main Toolbar")
        toolbar.setMovable(False)
        self.addToolBar(toolbar)

        # New file
        new_btn = QAction("New", self)
        new_btn.setToolTip("Create new input file (Ctrl+N)")
        new_btn.triggered.connect(self.new_file)
        toolbar.addAction(new_btn)

        # Open file
        open_btn = QAction("Open", self)
        open_btn.setToolTip("Open input file (Ctrl+O)")
        open_btn.triggered.connect(self.open_file)
        toolbar.addAction(open_btn)

        # Save file
        save_btn = QAction("Save", self)
        save_btn.setToolTip("Save input file (Ctrl+S)")
        save_btn.triggered.connect(self.save_file)
        toolbar.addAction(save_btn)

        toolbar.addSeparator()

        # Run SMOOTHIE - Create as QToolButton with custom styling
        run_tool_btn = QToolButton()
        run_tool_btn.setText("Run")
        run_tool_btn.setToolTip("Run SMOOTHIE calculation (Ctrl+R)")
        run_tool_btn.setObjectName("runButton")
        run_tool_btn.clicked.connect(self.run_smoothie)
        toolbar.addWidget(run_tool_btn)
        self.run_btn_widget = run_tool_btn

        # Stop SMOOTHIE - Create as QToolButton with custom styling
        stop_tool_btn = QToolButton()
        stop_tool_btn.setText("Stop")
        stop_tool_btn.setToolTip("Stop running calculation (Ctrl+.)")
        stop_tool_btn.setObjectName("stopButton")
        stop_tool_btn.setEnabled(False)
        stop_tool_btn.clicked.connect(self.stop_smoothie)
        toolbar.addWidget(stop_tool_btn)
        self.stop_btn_widget = stop_tool_btn

        # Keep QAction references for menu
        self.run_btn = QAction("Run", self)
        self.run_btn.setToolTip("Run SMOOTHIE calculation (Ctrl+R)")
        self.run_btn.triggered.connect(self.run_smoothie)

        self.stop_btn = QAction("Stop", self)
        self.stop_btn.setToolTip("Stop running calculation (Ctrl+.)")
        self.stop_btn.setEnabled(False)
        self.stop_btn.triggered.connect(self.stop_smoothie)

    def setup_connections(self):
        """Setup signal-slot connections"""
        self.smoothie_runner.output_ready.connect(self.log_widget.append_output)
        self.smoothie_runner.error_ready.connect(self.log_widget.append_error)
        self.smoothie_runner.finished.connect(self.on_calculation_finished)
        self.smoothie_runner.started.connect(self.on_calculation_started)

    def new_file(self):
        """Create a new input file"""
        reply = QMessageBox.question(
            self, "New File",
            "Clear all input fields?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )

        if reply == QMessageBox.Yes:
            self.input_panel.reset_all()
            self.current_file = None
            self.status_bar.showMessage("New file created")

    def open_file(self):
        """Open an existing input file"""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Open SMOOTHIE Input File",
            os.path.expanduser("~/Desktop/code/smoothie-public/smoothie/test"),
            "Input Files (*.in);;All Files (*)"
        )

        if file_path:
            try:
                self.input_panel.load_from_file(file_path)
                self.current_file = file_path
                self.status_bar.showMessage(f"Loaded: {os.path.basename(file_path)}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to open file:\n{str(e)}")

    def save_file(self):
        """Save the current input file"""
        if self.current_file:
            self.save_to_file(self.current_file)
        else:
            self.save_file_as()

    def save_file_as(self):
        """Save the input file with a new name"""
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save SMOOTHIE Input File",
            os.path.expanduser("~/Desktop/code/smoothie-public/smoothie/test"),
            "Input Files (*.in);;All Files (*)"
        )

        if file_path:
            self.save_to_file(file_path)
            self.current_file = file_path

    def save_to_file(self, file_path):
        """Save input data to a file"""
        try:
            content = self.input_panel.generate_input()
            with open(file_path, 'w') as f:
                f.write(content)
            self.status_bar.showMessage(f"Saved: {os.path.basename(file_path)}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save file:\n{str(e)}")

    def load_example(self):
        """Load an example input file"""
        example_path = os.path.expanduser(
            "~/Desktop/code/smoothie-public/smoothie/test/test.in"
        )

        if os.path.exists(example_path):
            try:
                self.input_panel.load_from_file(example_path)
                self.status_bar.showMessage("Example loaded successfully")
            except Exception as e:
                QMessageBox.warning(self, "Warning", f"Failed to load example:\n{str(e)}")
        else:
            QMessageBox.warning(self, "Warning", "Example file not found")

    def run_smoothie(self):
        """Run SMOOTHIE calculation"""
        # Ask user for working directory
        default_dir = self.working_directory if self.working_directory else os.path.expanduser(
            "~/Desktop/code/smoothie-public/smoothie/test"
        )

        work_dir = QFileDialog.getExistingDirectory(
            self,
            "Select Working Directory for Calculation",
            default_dir,
            QFileDialog.ShowDirsOnly | QFileDialog.DontResolveSymlinks
        )

        if not work_dir:
            # User cancelled
            return

        # Store the working directory
        self.working_directory = work_dir

        # Generate input file
        input_content = self.input_panel.generate_input()

        # Switch to log tab
        self.right_tabs.setCurrentWidget(self.log_widget)

        # Clear previous log
        self.log_widget.clear()
        self.log_widget.append_info("Starting SMOOTHIE calculation...")
        self.log_widget.append_info(f"Working directory: {work_dir}")

        # Save input file in the working directory
        input_file = os.path.join(work_dir, "smoothie_input.in")
        try:
            with open(input_file, 'w') as f:
                f.write(input_content)
            self.log_widget.append_info(f"Input file saved: {input_file}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to create input file:\n{str(e)}")
            return

        # Find SMOOTHIE executable
        smoothie_exe = os.path.expanduser(
            "~/Desktop/code/smoothie-public/smoothie/smoothie"
        )

        if not os.path.exists(smoothie_exe):
            QMessageBox.critical(
                self, "Error",
                f"SMOOTHIE executable not found at:\n{smoothie_exe}\n\n"
                "Please build SMOOTHIE first."
            )
            return

        # Pass working directory to runner
        self.is_running_smoothie = True  # Set flag before running SMOOTHIE
        self.smoothie_runner.run(smoothie_exe, input_file, work_dir)

    def stop_smoothie(self):
        """Stop the running SMOOTHIE calculation"""
        self.smoothie_runner.stop()
        self.log_widget.append_warning("Calculation stopped by user")

    def run_cm2lab(self):
        """Run CM2LAB for center-of-mass to lab frame conversion"""
        if not self.working_directory:
            QMessageBox.warning(
                self, "Warning",
                "Please run SMOOTHIE first to generate output files."
            )
            return

        # Generate cm2lab.in from SMOOTHIE input parameters
        try:
            cm2lab_input = self.generate_cm2lab_input()
            cm2lab_input_file = os.path.join(self.working_directory, "cm2lab.in")

            with open(cm2lab_input_file, 'w') as f:
                f.write(cm2lab_input)

            self.log_widget.append_info("Generated cm2lab.in")
            self.log_widget.append_info("Running CM2LAB...")

            # Find cm2lab executable
            cm2lab_exe = os.path.expanduser(
                "~/Desktop/code/smoothie-public/cm2lab/cm2lab"
            )

            if not os.path.exists(cm2lab_exe):
                QMessageBox.critical(
                    self, "Error",
                    f"CM2LAB executable not found at:\n{cm2lab_exe}\n\n"
                    "Please build CM2LAB first."
                )
                return

            # Run CM2LAB
            self.is_running_smoothie = False  # Set flag before running CM2LAB
            self.smoothie_runner.run(cm2lab_exe, cm2lab_input_file, self.working_directory)
            self.status_bar.showMessage("Running CM2LAB...")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to run CM2LAB:\n{str(e)}")
            self.log_widget.append_error(f"CM2LAB error: {str(e)}")

    def generate_cm2lab_input(self):
        """Generate cm2lab.in from current SMOOTHIE parameters"""
        # Get parameters from input panel
        m1 = self.input_panel.massp.value()  # projectile
        m2 = self.input_panel.masst.value()  # target
        m3 = self.input_panel.massb.value()  # detected
        m4 = m2 + self.input_panel.massx.value()  # recoil (target + undetected)

        elab = self.input_panel.elab.value()
        ecmfmin = self.input_panel.ecmbmin.value()
        ecmfmax = self.input_panel.ecmbmax.value()
        ecmfh = self.input_panel.ecmbh.value()

        thcmmin = self.input_panel.thmin.value()
        thcmmax = self.input_panel.thmax.value()
        thcminc = self.input_panel.thinc.value()

        lmax = self.input_panel.lxmax.value()

        # Generate cm2lab input file
        cm2lab_input = f"""NAMELIST
&cmsys  m1={m1:.1f} m2={m2:.1f} m3={m3:.1f} m4={m4:.1f} elab={elab:.1f} ecmfmin={ecmfmin:.1f}
          ecmfmax={ecmfmax:.1f} ecmfh={ecmfh:.1f} thcmmin={thcmmin:.1f} thcmmax={thcmmax:.1f} thcminc={thcminc:.1f} lmax={lmax} /
&labsys      elabmin=1 elabmax={int(elab+10)} elabh=2 thlabmin=1 thlabmax=180
          thlabinc=1 ptheta=2  pelab=2 /
"""

        return cm2lab_input

    def on_calculation_started(self):
        """Handle calculation start"""
        self.run_action.setEnabled(False)
        self.run_btn.setEnabled(False)
        self.run_btn_widget.setEnabled(False)
        self.stop_action.setEnabled(True)
        self.stop_btn.setEnabled(True)
        self.stop_btn_widget.setEnabled(True)
        self.status_bar.showMessage("Running SMOOTHIE...")

    def on_calculation_finished(self, exit_code):
        """Handle calculation completion"""
        self.run_action.setEnabled(True)
        self.run_btn.setEnabled(True)
        self.run_btn_widget.setEnabled(True)
        self.stop_action.setEnabled(False)
        self.stop_btn.setEnabled(False)
        self.stop_btn_widget.setEnabled(False)

        if exit_code == 0:
            if self.is_running_smoothie:
                # SMOOTHIE just finished successfully
                self.status_bar.showMessage("SMOOTHIE calculation completed successfully")
                self.log_widget.append_info("SMOOTHIE calculation completed successfully")

                # Set ecmb values for multi-group plots
                ecmbmin = self.input_panel.ecmbmin.value()
                ecmbmax = self.input_panel.ecmbmax.value()
                ecmbh = self.input_panel.ecmbh.value()
                self.plot_widget.set_ecmb_values(ecmbmin, ecmbmax, ecmbh)

                # Try to load and plot results from working directory
                if self.working_directory:
                    self.plot_widget.load_results(self.working_directory)
                else:
                    self.plot_widget.load_results()

                # Automatically run CM2LAB after successful SMOOTHIE run
                self.log_widget.append_info("Starting CM2LAB frame conversion automatically...")
                self.run_cm2lab()
            else:
                # CM2LAB just finished successfully
                self.status_bar.showMessage("CM2LAB completed successfully")
                self.log_widget.append_info("CM2LAB frame conversion completed successfully")

                # Set ecmb values for multi-group plots (in case needed for CM2LAB output)
                ecmbmin = self.input_panel.ecmbmin.value()
                ecmbmax = self.input_panel.ecmbmax.value()
                ecmbh = self.input_panel.ecmbh.value()
                self.plot_widget.set_ecmb_values(ecmbmin, ecmbmax, ecmbh)

                # Reload results to include CM2LAB output
                if self.working_directory:
                    self.plot_widget.load_results(self.working_directory)
        else:
            self.status_bar.showMessage(f"Calculation failed (exit code: {exit_code})")
            self.log_widget.append_error(f"Calculation failed with exit code: {exit_code}")

    def change_theme(self, theme):
        """Change the application theme"""
        apply_modern_style(self, theme)
        self.status_bar.showMessage(f"Theme changed to {theme}")

    def show_about(self):
        """Show about dialog"""
        QMessageBox.about(
            self,
            "About SMOOTHIE",
            "<h2>SMOOTHIE</h2>"
            "<p><b>S</b>cattering <b>M</b>odel of <b>O</b>ptical <b>O</b>perator "
            "<b>T</b>heory for <b>I</b>chimura-Austern-Vincent <b>E</b>quations</p>"
            "<p><b>IAV nonelastic breakup calculations</b></p>"
            "<p>Modern GUI for quantum scattering calculations</p>"
            "<p>Version 1.1</p>"
            "<p>Built with PySide6</p>"
        )

    def show_docs(self):
        """Show documentation"""
        QMessageBox.information(
            self,
            "Documentation",
            "For documentation, please visit:\n\n"
            "https://smoothie.fewbody.com\n\n"
            "Or check the README files in the source directory."
        )
