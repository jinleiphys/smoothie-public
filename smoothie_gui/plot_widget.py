"""
Plot widget for displaying SMOOTHIE results using matplotlib
"""

from PySide6.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QComboBox, QLabel
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.interpolate import griddata


class PlotWidget(QWidget):
    """Widget for plotting SMOOTHIE results"""

    def __init__(self):
        super().__init__()
        self.current_data = {}
        self.working_directory = None  # Store working directory for refresh
        self.ecmb_values = []  # Store ecmb energy values for labels
        self.init_ui()

    def init_ui(self):
        """Initialize the plot widget"""
        layout = QVBoxLayout(self)

        # Control panel
        control_layout = QHBoxLayout()

        self.plot_type = QComboBox()
        self.plot_type.addItems([
            "Double X-Section Angular (fort.21)",
            "X-Section Energy (fort.22)",
            "X-Section Alpha_xA (fort.20)",
            "X-Section La (fort.23)",
            "X-Section Lb (fort.24)",
            "X-Section Energy Distribution (fort.912)",
            "X-Section Angular Distribution (fort.913)",
            "Double X-Section 3D Lab (fort.917)"
        ])
        self.plot_type.currentIndexChanged.connect(self.update_plot)
        control_layout.addWidget(QLabel("Plot Type:"))
        control_layout.addWidget(self.plot_type)

        control_layout.addStretch()

        refresh_btn = QPushButton("Refresh")
        refresh_btn.clicked.connect(self.load_results)
        control_layout.addWidget(refresh_btn)

        clear_btn = QPushButton("Clear")
        clear_btn.clicked.connect(self.clear_plot)
        control_layout.addWidget(clear_btn)

        layout.addLayout(control_layout)

        # Matplotlib figure
        self.figure = Figure(figsize=(8, 6), dpi=100)
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)

        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)

        # Initial plot
        self.clear_plot()

    def set_ecmb_values(self, ecmbmin, ecmbmax, ecmbh):
        """Set the ecmb energy values for multi-group plots"""
        self.ecmb_values = []
        ecmb = ecmbmin
        while ecmb <= ecmbmax + 1e-6:  # Small epsilon for floating point comparison
            self.ecmb_values.append(ecmb)
            ecmb += ecmbh

    def load_results(self, working_dir=None):
        """Load results from SMOOTHIE output files"""
        try:
            # Store working directory for future refreshes
            if working_dir is not None:
                self.working_directory = working_dir

            # Use stored working directory if available
            if self.working_directory is None:
                self.working_directory = os.getcwd()

            # Try to read common output files
            # fort.20 - cross sections
            # fort.21 - angular distributions
            # fort.22 - energy spectra

            # fort.20, 21, 23, 24 have multiple groups separated by &
            fort20 = os.path.join(self.working_directory, 'fort.20')
            if os.path.exists(fort20):
                data = self.read_fort_file_multigroup(fort20)
                self.current_data['fort20'] = data

            fort21 = os.path.join(self.working_directory, 'fort.21')
            if os.path.exists(fort21):
                data = self.read_fort_file_multigroup(fort21)
                self.current_data['fort21'] = data

            fort22 = os.path.join(self.working_directory, 'fort.22')
            if os.path.exists(fort22):
                data = self.read_fort_file(fort22)
                self.current_data['fort22'] = data

            fort23 = os.path.join(self.working_directory, 'fort.23')
            if os.path.exists(fort23):
                data = self.read_fort_file_multigroup(fort23)
                self.current_data['fort23'] = data

            fort24 = os.path.join(self.working_directory, 'fort.24')
            if os.path.exists(fort24):
                data = self.read_fort_file_multigroup(fort24)
                self.current_data['fort24'] = data

            fort912 = os.path.join(self.working_directory, 'fort.912')
            if os.path.exists(fort912):
                data = self.read_fort_file(fort912)
                self.current_data['fort912'] = data

            fort913 = os.path.join(self.working_directory, 'fort.913')
            if os.path.exists(fort913):
                data = self.read_fort_file(fort913)
                self.current_data['fort913'] = data

            fort917 = os.path.join(self.working_directory, 'fort.917')
            if os.path.exists(fort917):
                data = self.read_fort_file(fort917)
                self.current_data['fort917'] = data

            self.update_plot()

        except Exception as e:
            print(f"Error loading results: {e}")

    def read_fort_file(self, filename):
        """Read a FORTRAN output file"""
        try:
            data = np.loadtxt(filename)
            return data
        except:
            # Try reading with more flexible parsing
            with open(filename, 'r') as f:
                lines = f.readlines()

            data_lines = []
            for line in lines:
                line = line.strip()
                if line and not line.startswith('#'):
                    try:
                        values = [float(x) for x in line.split()]
                        data_lines.append(values)
                    except:
                        continue

            if data_lines:
                return np.array(data_lines)
            return None

    def read_fort_file_multigroup(self, filename):
        """Read a FORTRAN output file with multiple groups separated by &"""
        try:
            with open(filename, 'r') as f:
                lines = f.readlines()

            groups = []
            current_group = []

            for line in lines:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                # Check if this line is a separator
                if line.startswith('&'):
                    if current_group:
                        groups.append(np.array(current_group))
                        current_group = []
                else:
                    try:
                        values = [float(x) for x in line.split()]
                        current_group.append(values)
                    except:
                        continue

            # Add the last group
            if current_group:
                groups.append(np.array(current_group))

            return groups if groups else None
        except Exception as e:
            print(f"Error reading file {filename}: {e}")
            return None

    def update_plot(self):
        """Update the plot based on selected type"""
        self.figure.clear()
        ax = self.figure.add_subplot(111)

        plot_idx = self.plot_type.currentIndex()

        if plot_idx == 0:  # Double X-Section Angular (fort.21)
            self.plot_fort21(ax)
        elif plot_idx == 1:  # X-Section Energy (fort.22)
            self.plot_fort22(ax)
        elif plot_idx == 2:  # X-Section Alpha_xA (fort.20)
            self.plot_fort20(ax)
        elif plot_idx == 3:  # X-Section La (fort.23)
            self.plot_fort23(ax)
        elif plot_idx == 4:  # X-Section Lb (fort.24)
            self.plot_fort24(ax)
        elif plot_idx == 5:  # X-Section Energy Distribution (fort.912)
            self.plot_fort912(ax)
        elif plot_idx == 6:  # X-Section Angular Distribution (fort.913)
            self.plot_fort913(ax)
        elif plot_idx == 7:  # Double X-Section 3D Lab (fort.917)
            self.plot_fort917(ax)

        self.canvas.draw()

    def plot_fort20(self, ax):
        """Plot fort.20: x-section alpha_xA distribution (multiple groups)"""
        if 'fort20' in self.current_data:
            groups = self.current_data['fort20']
            if groups is not None and isinstance(groups, list) and len(groups) > 0:
                # Color cycle for multiple groups
                colors = plt.cm.tab10(np.linspace(0, 1, max(len(groups), 10)))

                for i, data in enumerate(groups):
                    if len(data.shape) >= 2 and data.shape[1] >= 2:
                        color = colors[i % len(colors)]
                        # Use ecmb value for label if available
                        if i < len(self.ecmb_values):
                            label = f'Ecmb = {self.ecmb_values[i]:.1f} MeV'
                        else:
                            label = f'Group {i+1}'
                        ax.plot(data[:, 0], data[:, 1], '-', color=color, linewidth=2,
                               marker='o', markersize=4, label=label)

                ax.set_xlabel('Alpha_xA', fontsize=12)
                ax.set_ylabel('Cross Section (mb)', fontsize=12)
                ax.set_title('X-Section Alpha_xA Distribution (fort.20)', fontsize=14, fontweight='bold')
                ax.grid(True, alpha=0.3)
                if len(groups) > 1:
                    ax.legend()
            else:
                self._show_no_data(ax, 'fort.20')
        else:
            self._show_no_data(ax, 'fort.20')

    def plot_fort21(self, ax):
        """Plot fort.21: double x-section angular distribution (multiple groups)"""
        if 'fort21' in self.current_data:
            groups = self.current_data['fort21']
            if groups is not None and isinstance(groups, list) and len(groups) > 0:
                # Color cycle for multiple groups
                colors = plt.cm.tab10(np.linspace(0, 1, max(len(groups), 10)))

                for i, data in enumerate(groups):
                    if len(data.shape) >= 2 and data.shape[1] >= 2:
                        color = colors[i % len(colors)]
                        # Use ecmb value for label if available
                        if i < len(self.ecmb_values):
                            label = f'Ecmb = {self.ecmb_values[i]:.1f} MeV'
                        else:
                            label = f'Group {i+1}'
                        ax.plot(data[:, 0], data[:, 1], '-', color=color, linewidth=2,
                               marker='s', markersize=4, label=label)

                ax.set_xlabel('Angle (degrees)', fontsize=12)
                ax.set_ylabel('d²σ/dΩdE (mb/sr/MeV)', fontsize=12)
                ax.set_title('Double X-Section Angular Distribution (fort.21)', fontsize=14, fontweight='bold')
                ax.grid(True, alpha=0.3)
                ax.set_yscale('log')
                if len(groups) > 1:
                    ax.legend()
            else:
                self._show_no_data(ax, 'fort.21')
        else:
            self._show_no_data(ax, 'fort.21')

    def plot_fort22(self, ax):
        """Plot fort.22: x-section energy distribution"""
        if 'fort22' in self.current_data:
            data = self.current_data['fort22']
            if data is not None and len(data.shape) >= 2 and data.shape[1] >= 2:
                ax.plot(data[:, 0], data[:, 1], 'k-', linewidth=2, marker='^', markersize=4)
                ax.set_xlabel('Energy (MeV)', fontsize=12)
                ax.set_ylabel('dσ/dE (mb/MeV)', fontsize=12)
                ax.set_title('X-Section Energy Distribution (fort.22)', fontsize=14, fontweight='bold')
                ax.grid(True, alpha=0.3)
            else:
                self._show_no_data(ax, 'fort.22')
        else:
            self._show_no_data(ax, 'fort.22')

    def plot_fort23(self, ax):
        """Plot fort.23: x-section la distribution (multiple groups)"""
        if 'fort23' in self.current_data:
            groups = self.current_data['fort23']
            if groups is not None and isinstance(groups, list) and len(groups) > 0:
                # Color cycle for multiple groups
                colors = plt.cm.tab10(np.linspace(0, 1, max(len(groups), 10)))

                for i, data in enumerate(groups):
                    if len(data.shape) >= 2 and data.shape[1] >= 2:
                        color = colors[i % len(colors)]
                        # Use ecmb value for label if available
                        if i < len(self.ecmb_values):
                            label = f'Ecmb = {self.ecmb_values[i]:.1f} MeV'
                        else:
                            label = f'Group {i+1}'
                        ax.plot(data[:, 0], data[:, 1], '-', color=color, linewidth=2,
                               marker='o', markersize=6, label=label)

                ax.set_xlabel('La (angular momentum)', fontsize=12)
                ax.set_ylabel('Cross Section (mb)', fontsize=12)
                ax.set_title('X-Section La Distribution (fort.23)', fontsize=14, fontweight='bold')
                ax.grid(True, alpha=0.3)
                if len(groups) > 1:
                    ax.legend()
            else:
                self._show_no_data(ax, 'fort.23')
        else:
            self._show_no_data(ax, 'fort.23')

    def plot_fort24(self, ax):
        """Plot fort.24: x-section lb distribution (multiple groups)"""
        if 'fort24' in self.current_data:
            groups = self.current_data['fort24']
            if groups is not None and isinstance(groups, list) and len(groups) > 0:
                # Color cycle for multiple groups
                colors = plt.cm.tab10(np.linspace(0, 1, max(len(groups), 10)))

                for i, data in enumerate(groups):
                    if len(data.shape) >= 2 and data.shape[1] >= 2:
                        color = colors[i % len(colors)]
                        # Use ecmb value for label if available
                        if i < len(self.ecmb_values):
                            label = f'Ecmb = {self.ecmb_values[i]:.1f} MeV'
                        else:
                            label = f'Group {i+1}'
                        ax.plot(data[:, 0], data[:, 1], '-', color=color, linewidth=2,
                               marker='s', markersize=6, label=label)

                ax.set_xlabel('Lb (angular momentum)', fontsize=12)
                ax.set_ylabel('Cross Section (mb)', fontsize=12)
                ax.set_title('X-Section Lb Distribution (fort.24)', fontsize=14, fontweight='bold')
                ax.grid(True, alpha=0.3)
                if len(groups) > 1:
                    ax.legend()
            else:
                self._show_no_data(ax, 'fort.24')
        else:
            self._show_no_data(ax, 'fort.24')

    def plot_fort912(self, ax):
        """Plot fort.912: cross section energy distribution"""
        if 'fort912' in self.current_data:
            data = self.current_data['fort912']
            if data is not None and len(data.shape) >= 2 and data.shape[1] >= 2:
                ax.plot(data[:, 0], data[:, 1], 'k-', linewidth=2, marker='o', markersize=4)
                ax.set_xlabel('Energy (MeV)', fontsize=12)
                ax.set_ylabel('Cross Section (mb/MeV)', fontsize=12)
                ax.set_title('Cross Section Energy Distribution (fort.912)', fontsize=14, fontweight='bold')
                ax.grid(True, alpha=0.3)
            else:
                self._show_no_data(ax, 'fort.912')
        else:
            self._show_no_data(ax, 'fort.912')

    def plot_fort913(self, ax):
        """Plot fort.913: cross section angular distribution"""
        if 'fort913' in self.current_data:
            data = self.current_data['fort913']
            if data is not None and len(data.shape) >= 2 and data.shape[1] >= 2:
                ax.plot(data[:, 0], data[:, 1], 'k-', linewidth=2, marker='s', markersize=4)
                ax.set_xlabel('Angle (degrees)', fontsize=12)
                ax.set_ylabel('Cross Section (mb/sr)', fontsize=12)
                ax.set_title('Cross Section Angular Distribution (fort.913)', fontsize=14, fontweight='bold')
                ax.grid(True, alpha=0.3)
            else:
                self._show_no_data(ax, 'fort.913')
        else:
            self._show_no_data(ax, 'fort.913')

    def plot_fort917(self, ax):
        """Plot fort.917: double cross sections 3-D plot in lab (contour)"""
        if 'fort917' in self.current_data:
            data = self.current_data['fort917']
            if data is not None and len(data.shape) >= 2 and data.shape[1] >= 3:
                # Assuming data format: angle, energy, cross_section
                try:
                    # Extract data points
                    angles = data[:, 0]
                    energies = data[:, 1]
                    cross_sections = data[:, 2]

                    # Create a finer grid for interpolation
                    angles_unique = np.unique(angles)
                    energies_unique = np.unique(energies)

                    # Create a much denser grid (3x resolution)
                    angles_interp = np.linspace(angles_unique.min(), angles_unique.max(), len(angles_unique) * 3)
                    energies_interp = np.linspace(energies_unique.min(), energies_unique.max(), len(energies_unique) * 3)
                    X_interp, Y_interp = np.meshgrid(angles_interp, energies_interp)

                    # Interpolate data onto the finer grid
                    points = np.column_stack((angles, energies))
                    Z_interp = griddata(points, cross_sections, (X_interp, Y_interp), method='cubic')

                    # Create smooth contour plot with more levels
                    contour = ax.contourf(X_interp, Y_interp, Z_interp, levels=50, cmap='viridis', antialiased=True)
                    self.figure.colorbar(contour, ax=ax, label='d²σ/dΩdE (mb/sr/MeV)')

                    ax.set_xlabel('Angle (degrees)', fontsize=12)
                    ax.set_ylabel('Energy (MeV)', fontsize=12)
                    ax.set_title('Double X-Section 3D Lab Frame (fort.917)', fontsize=14, fontweight='bold')
                except Exception as e:
                    print(f"Error creating contour plot for fort.917: {e}")
                    self._show_no_data(ax, 'fort.917')
            else:
                self._show_no_data(ax, 'fort.917')
        else:
            self._show_no_data(ax, 'fort.917')

    def _show_no_data(self, ax, filename):
        """Display no data message"""
        ax.text(0.5, 0.5, f'No data available for {filename}\n\nRun SMOOTHIE to generate results',
               ha='center', va='center', transform=ax.transAxes,
               fontsize=12, color='gray')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')

    def clear_plot(self):
        """Clear the plot"""
        self.current_data = {}
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.text(0.5, 0.5, 'No data to display\n\nRun a calculation to see results',
               ha='center', va='center', transform=ax.transAxes,
               fontsize=14, color='gray')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        self.canvas.draw()
