"""
Input panel with tabs for different parameter sections
"""

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QTabWidget, QFormLayout, QLineEdit,
    QSpinBox, QDoubleSpinBox, QComboBox, QGroupBox, QScrollArea, QLabel
)
from PySide6.QtCore import Qt
import re


class InputPanel(QWidget):
    """Panel containing all input forms organized in tabs"""

    def __init__(self):
        super().__init__()
        self.init_ui()

    def init_ui(self):
        """Initialize the UI with tabs for different sections"""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(10, 10, 10, 10)

        # Create tab widget
        self.tabs = QTabWidget()
        self.tabs.setDocumentMode(True)

        # Create tabs for each section
        self.global_tab = self.create_global_tab()
        self.system_tab = self.create_system_tab()
        self.outgoing_tab = self.create_outgoing_tab()
        self.potential_tab = self.create_potential_tab()

        self.tabs.addTab(self.global_tab, "Global")
        self.tabs.addTab(self.system_tab, "System")
        self.tabs.addTab(self.outgoing_tab, "Outgoing")
        self.tabs.addTab(self.potential_tab, "Potentials")

        layout.addWidget(self.tabs)

    def create_global_tab(self):
        """Create the GLOBAL parameters tab"""
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.NoFrame)

        widget = QWidget()
        layout = QFormLayout(widget)

        # Radial parameters
        self.hcm = QDoubleSpinBox()
        self.hcm.setRange(0.001, 1.0)
        self.hcm.setSingleStep(0.001)
        self.hcm.setDecimals(3)
        self.hcm.setValue(0.05)
        layout.addRow("Radial Step Size (fm):", self.hcm)

        self.rmax = QDoubleSpinBox()
        self.rmax.setRange(1.0, 1000.0)
        self.rmax.setSingleStep(1.0)
        self.rmax.setValue(50.0)
        layout.addRow("Maximum Radius (fm):", self.rmax)

        # Angular momentum parameters
        self.lmin = QSpinBox()
        self.lmin.setRange(0, 100)
        self.lmin.setValue(0)
        layout.addRow("Minimum L (λₐ, λᵦ):", self.lmin)

        self.lmax = QSpinBox()
        self.lmax.setRange(1, 100)
        self.lmax.setValue(25)
        layout.addRow("Maximum L (λₐ, λᵦ):", self.lmax)

        self.jtmin = QSpinBox()
        self.jtmin.setRange(0, 100)
        self.jtmin.setValue(0)
        layout.addRow("Minimum Total J:", self.jtmin)

        self.jtmax = QSpinBox()
        self.jtmax.setRange(1, 100)
        self.jtmax.setValue(25)
        layout.addRow("Maximum Total J:", self.jtmax)

        self.lxmin = QSpinBox()
        self.lxmin.setRange(0, 100)
        self.lxmin.setValue(0)
        layout.addRow("Minimum lₓ:", self.lxmin)

        self.lxmax = QSpinBox()
        self.lxmax.setRange(1, 100)
        self.lxmax.setValue(12)
        layout.addRow("Maximum lₓ:", self.lxmax)

        # Energy and method
        self.elab = QDoubleSpinBox()
        self.elab.setRange(0.1, 1000.0)
        self.elab.setSingleStep(0.1)
        self.elab.setValue(25.5)
        layout.addRow("Incident Energy (MeV):", self.elab)

        self.dwba = QComboBox()
        self.dwba.addItems([
            "1: No spins, rᵦₓ variable",
            "2: No spins, rᵦ variable",
            "3: With spins, rᵦ variable",
            "4: Enhanced Lagrange mesh",
            "5: Lagrange mesh, rᵦ variable"
        ])
        layout.addRow("DWBA Method:", self.dwba)

        # Angular range
        self.thmin = QDoubleSpinBox()
        self.thmin.setRange(0.0, 180.0)
        self.thmin.setValue(0.0)
        layout.addRow("Minimum Angle (degrees):", self.thmin)

        self.thmax = QDoubleSpinBox()
        self.thmax.setRange(0.0, 180.0)
        self.thmax.setValue(180.0)
        layout.addRow("Maximum Angle (degrees):", self.thmax)

        self.thinc = QDoubleSpinBox()
        self.thinc.setRange(0.1, 10.0)
        self.thinc.setValue(1.0)
        layout.addRow("Angular Increment (degrees):", self.thinc)

        # Quadrature points
        self.nx = QSpinBox()
        self.nx.setRange(1, 200)
        self.nx.setValue(34)
        layout.addRow("Gaussian Quadrature Points (cos θ):", self.nx)

        self.nr = QSpinBox()
        self.nr.setRange(1, 1000)
        self.nr.setValue(100)
        layout.addRow("Radial Quadrature Points:", self.nr)

        # Output options
        self.printf = QComboBox()
        self.printf.addItems(["False", "True"])
        layout.addRow("Enable Detailed Output:", self.printf)

        scroll.setWidget(widget)
        return scroll

    def create_system_tab(self):
        """Create the SYSTEM parameters tab"""
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.NoFrame)

        widget = QWidget()
        layout = QFormLayout(widget)

        # Projectile
        group_p = QGroupBox("Projectile (a)")
        layout_p = QFormLayout()
        self.namep = QLineEdit("d")
        layout_p.addRow("Name:", self.namep)
        self.massp = QDoubleSpinBox()
        self.massp.setRange(0.001, 300.0)
        self.massp.setDecimals(4)
        self.massp.setValue(2.0)
        layout_p.addRow("Mass (amu):", self.massp)
        self.zp = QDoubleSpinBox()
        self.zp.setRange(0.0, 100.0)
        self.zp.setValue(1.0)
        layout_p.addRow("Charge:", self.zp)
        self.jp = QDoubleSpinBox()
        self.jp.setRange(0.0, 20.0)
        self.jp.setSingleStep(0.5)
        self.jp.setValue(0.0)
        layout_p.addRow("Spin:", self.jp)
        group_p.setLayout(layout_p)
        layout.addRow(group_p)

        # Target
        group_t = QGroupBox("Target (A)")
        layout_t = QFormLayout()
        self.namet = QLineEdit("93Nb")
        layout_t.addRow("Name:", self.namet)
        self.masst = QDoubleSpinBox()
        self.masst.setRange(0.001, 300.0)
        self.masst.setDecimals(4)
        self.masst.setValue(93.0)
        layout_t.addRow("Mass (amu):", self.masst)
        self.zt = QDoubleSpinBox()
        self.zt.setRange(0.0, 100.0)
        self.zt.setValue(41.0)
        layout_t.addRow("Charge:", self.zt)
        self.jt = QDoubleSpinBox()
        self.jt.setRange(0.0, 20.0)
        self.jt.setSingleStep(0.5)
        self.jt.setValue(0.0)
        layout_t.addRow("Spin:", self.jt)
        group_t.setLayout(layout_t)
        layout.addRow(group_t)

        # Detected particle
        group_b = QGroupBox("Detected Particle (b)")
        layout_b = QFormLayout()
        self.nameb = QLineEdit("p")
        layout_b.addRow("Name:", self.nameb)
        self.massb = QDoubleSpinBox()
        self.massb.setRange(0.001, 300.0)
        self.massb.setDecimals(4)
        self.massb.setValue(1.0078)
        layout_b.addRow("Mass (amu):", self.massb)
        self.zb = QDoubleSpinBox()
        self.zb.setRange(0.0, 100.0)
        self.zb.setValue(1.0)
        layout_b.addRow("Charge:", self.zb)
        self.jb = QDoubleSpinBox()
        self.jb.setRange(0.0, 20.0)
        self.jb.setSingleStep(0.5)
        self.jb.setValue(0.0)
        layout_b.addRow("Spin:", self.jb)
        group_b.setLayout(layout_b)
        layout.addRow(group_b)

        # Undetected particle
        group_x = QGroupBox("Undetected Particle (x)")
        layout_x = QFormLayout()
        self.namex = QLineEdit("n")
        layout_x.addRow("Name:", self.namex)
        self.massx = QDoubleSpinBox()
        self.massx.setRange(0.001, 300.0)
        self.massx.setDecimals(4)
        self.massx.setValue(1.0087)
        layout_x.addRow("Mass (amu):", self.massx)
        self.zx = QDoubleSpinBox()
        self.zx.setRange(0.0, 100.0)
        self.zx.setValue(0.0)
        layout_x.addRow("Charge:", self.zx)
        self.jx = QDoubleSpinBox()
        self.jx.setRange(0.0, 20.0)
        self.jx.setSingleStep(0.5)
        self.jx.setValue(0.0)
        layout_x.addRow("Spin:", self.jx)
        group_x.setLayout(layout_x)
        layout.addRow(group_x)

        # System properties
        group_sys = QGroupBox("System Properties")
        layout_sys = QFormLayout()
        self.sbx = QDoubleSpinBox()
        self.sbx.setRange(0.0, 20.0)
        self.sbx.setSingleStep(0.5)
        self.sbx.setValue(0.0)
        layout_sys.addRow("Total Spin Coupling (b+x):", self.sbx)
        self.lbx = QSpinBox()
        self.lbx.setRange(0, 20)
        self.lbx.setValue(0)
        layout_sys.addRow("Orbital L (b+x):", self.lbx)
        self.nodes = QSpinBox()
        self.nodes.setRange(1, 10)
        self.nodes.setValue(1)
        layout_sys.addRow("Number of Nodes:", self.nodes)
        self.be = QDoubleSpinBox()
        self.be.setRange(0.001, 100.0)
        self.be.setDecimals(3)
        self.be.setValue(2.224)
        layout_sys.addRow("Binding Energy (MeV):", self.be)
        group_sys.setLayout(layout_sys)
        layout.addRow(group_sys)

        scroll.setWidget(widget)
        return scroll

    def create_outgoing_tab(self):
        """Create the OUTGOING parameters tab"""
        widget = QWidget()
        layout = QFormLayout(widget)

        self.ecmbmin = QDoubleSpinBox()
        self.ecmbmin.setRange(0.1, 1000.0)
        self.ecmbmin.setValue(2.0)
        layout.addRow("Minimum Energy (MeV):", self.ecmbmin)

        self.ecmbmax = QDoubleSpinBox()
        self.ecmbmax.setRange(0.1, 1000.0)
        self.ecmbmax.setValue(30.0)
        layout.addRow("Maximum Energy (MeV):", self.ecmbmax)

        self.ecmbh = QDoubleSpinBox()
        self.ecmbh.setRange(0.1, 100.0)
        self.ecmbh.setValue(1.0)
        layout.addRow("Energy Step (MeV):", self.ecmbh)

        return widget

    def create_potential_tab(self):
        """Create the POTENTIAL parameters tab"""
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.NoFrame)

        widget = QWidget()
        main_layout = QVBoxLayout(widget)

        # Create tabs for each potential
        pot_tabs = QTabWidget()

        self.pot1 = self.create_potential_group("a", "Projectile-Target")
        self.pot2 = self.create_potential_group("b", "Detected-Recoil")
        self.pot3 = self.create_potential_group("x", "Undetected-Target")
        self.pot4 = self.create_potential_group("p", "Detected-Undetected")
        self.pot5 = self.create_potential_group("t", "Detected-Target")

        pot_tabs.addTab(self.pot1['widget'], f"Pot 1 (a)")
        pot_tabs.addTab(self.pot2['widget'], f"Pot 2 (b)")
        pot_tabs.addTab(self.pot3['widget'], f"Pot 3 (x)")
        pot_tabs.addTab(self.pot4['widget'], f"Pot 4 (p)")
        pot_tabs.addTab(self.pot5['widget'], f"Pot 5 (t)")

        main_layout.addWidget(pot_tabs)

        scroll.setWidget(widget)
        return scroll

    def create_potential_group(self, kp1, description):
        """Create a potential parameter group"""
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.NoFrame)

        widget = QWidget()
        layout = QFormLayout(widget)

        pot_dict = {'kp1': kp1}

        # Potential type
        ptype = QComboBox()
        ptype.addItems([
            "1: Woods-Saxon",
            "2: Gaussian",
            "3: KD02 Global",
            "4: CH89 Global",
            "10: Global Deuteron",
            "41: Read from fort.41",
            "42: Read from fort.42",
            "43: Read from fort.43",
            "44: Read from fort.44"
        ])
        if kp1 in ['b', 'x', 't']:
            ptype.setCurrentIndex(3)  # CH89 by default
        elif kp1 == 'p':
            ptype.setCurrentIndex(1)  # Gaussian by default
        layout.addRow("Potential Model:", ptype)
        pot_dict['ptype'] = ptype

        # Mass parameters
        a1 = QDoubleSpinBox()
        a1.setRange(0.0, 300.0)
        a1.setDecimals(3)
        a1.setValue(0.0 if kp1 != 'p' else 1.0)
        layout.addRow("Mass Parameter 1 (amu):", a1)
        pot_dict['a1'] = a1

        a2 = QDoubleSpinBox()
        a2.setRange(0.0, 300.0)
        a2.setDecimals(3)
        a2.setValue(93.0 if kp1 in ['a', 'x', 't'] else (94.0 if kp1 == 'b' else 1.0))
        layout.addRow("Mass Parameter 2 (amu):", a2)
        pot_dict['a2'] = a2

        rc = QDoubleSpinBox()
        rc.setRange(0.0, 10.0)
        rc.setDecimals(2)
        rc.setValue(1.3 if kp1 == 'a' else (1.5 if kp1 == 'p' else 0.0))
        layout.addRow("Coulomb Radius (fm):", rc)
        pot_dict['rc'] = rc

        # Woods-Saxon parameters
        uv = QDoubleSpinBox()
        uv.setRange(-200.0, 200.0)
        uv.setDecimals(2)
        uv.setValue(77.3 if kp1 == 'a' else (72.15 if kp1 == 'p' else 0.0))
        layout.addRow("Real Volume Depth (MeV):", uv)
        pot_dict['uv'] = uv

        av = QDoubleSpinBox()
        av.setRange(0.0, 10.0)
        av.setDecimals(3)
        av.setValue(0.77 if kp1 == 'a' else (1.484 if kp1 == 'p' else 0.0))
        layout.addRow("Real Volume Diffuseness (fm):", av)
        pot_dict['av'] = av

        rv = QDoubleSpinBox()
        rv.setRange(0.0, 10.0)
        rv.setDecimals(2)
        rv.setValue(1.15 if kp1 == 'a' else 0.0)
        layout.addRow("Real Volume Radius (fm):", rv)
        pot_dict['rv'] = rv

        uw = QDoubleSpinBox()
        uw.setRange(-200.0, 200.0)
        uw.setDecimals(2)
        uw.setValue(6.1 if kp1 == 'a' else 0.0)
        layout.addRow("Imaginary Volume Depth (MeV):", uw)
        pot_dict['uw'] = uw

        aw = QDoubleSpinBox()
        aw.setRange(0.0, 10.0)
        aw.setDecimals(3)
        aw.setValue(0.47 if kp1 == 'a' else 0.0)
        layout.addRow("Imaginary Volume Diffuseness (fm):", aw)
        pot_dict['aw'] = aw

        rw = QDoubleSpinBox()
        rw.setRange(0.0, 10.0)
        rw.setDecimals(2)
        rw.setValue(1.33 if kp1 == 'a' else 0.0)
        layout.addRow("Imaginary Volume Radius (fm):", rw)
        pot_dict['rw'] = rw

        wd = QDoubleSpinBox()
        wd.setRange(-200.0, 200.0)
        wd.setDecimals(2)
        wd.setValue(8.4 if kp1 == 'a' else 0.0)
        layout.addRow("Surface Absorption Depth (MeV):", wd)
        pot_dict['wd'] = wd

        awd = QDoubleSpinBox()
        awd.setRange(0.0, 10.0)
        awd.setDecimals(3)
        awd.setValue(0.77 if kp1 == 'a' else 0.0)
        layout.addRow("Surface Absorption Diffuseness (fm):", awd)
        pot_dict['awd'] = awd

        rwd = QDoubleSpinBox()
        rwd.setRange(0.0, 10.0)
        rwd.setDecimals(2)
        rwd.setValue(1.37 if kp1 == 'a' else 0.0)
        layout.addRow("Surface Absorption Radius (fm):", rwd)
        pot_dict['rwd'] = rwd

        pot_dict['widget'] = scroll
        scroll.setWidget(widget)

        return pot_dict

    def generate_input(self):
        """Generate SMOOTHIE input file content"""
        lines = ["NAMELIST"]

        # GLOBAL section
        lines.append("&GLOBAL      hcm={:.3f}  lmax={}  elab={:.1f} thmin={:.1f} thmax={:.1f}  printf={} dwba={}".format(
            self.hcm.value(),
            self.lmax.value(),
            self.elab.value(),
            self.thmin.value(),
            self.thmax.value(),
            't' if self.printf.currentIndex() == 1 else 'f',
            self.dwba.currentIndex() + 1
        ))
        lines.append("             thinc={:.1f}   nx={}  rmax={:.1f}   nr={}  lxmax={}  /".format(
            self.thinc.value(),
            self.nx.value(),
            self.rmax.value(),
            self.nr.value(),
            self.lxmax.value()
        ))
        lines.append("\n\n")

        # SYSTEM section
        lines.append("&SYSTEM     namep='{}'     massp={:.1f}       zp={:.1f}    jp={:.1f} sbx={:.1f}".format(
            self.namep.text(),
            self.massp.value(),
            self.zp.value(),
            self.jp.value(),
            self.sbx.value()
        ))
        lines.append("            namet='{}'  masst={:.1f}     zt={:.1f}   jt={:.1f}  be={:.3f}".format(
            self.namet.text(),
            self.masst.value(),
            self.zt.value(),
            self.jt.value(),
            self.be.value()
        ))
        lines.append("            nameb='{}'     massb={:.4f}        zb={:.1f}    jb={:.1f}".format(
            self.nameb.text(),
            self.massb.value(),
            self.zb.value(),
            self.jb.value()
        ))
        lines.append("            namex='{}'     massx={:.4f}   zx={:.1f}    jx={:.1f}  lbx={}    nodes={}  /".format(
            self.namex.text(),
            self.massx.value(),
            self.zx.value(),
            self.jx.value(),
            self.lbx.value(),
            self.nodes.value()
        ))
        lines.append("\n\n")

        # OUTGOING section
        lines.append("&OUTGOING   ecmbmin={:.1f} ecmbmax={:.1f} ecmbh={:.1f}  /".format(
            self.ecmbmin.value(),
            self.ecmbmax.value(),
            self.ecmbh.value()
        ))
        lines.append("\n\n")
        lines.append("&OUTGOING /\n\n\n")

        # POTENTIAL sections
        for pot in [self.pot1, self.pot2, self.pot3, self.pot4, self.pot5]:
            lines.append(self.generate_potential_section(pot))

        lines.append("&POTENTIAL /\n\n\n\n\n")

        return '\n'.join(lines)

    def generate_potential_section(self, pot):
        """Generate a POTENTIAL namelist section"""
        ptype = pot['ptype'].currentIndex() + 1
        if ptype > 5:
            ptype = 41 + (ptype - 6)

        lines = ["&POTENTIAL  kp1='{}'  ptype={} a1={:.0f} a2={:.0f}".format(
            pot['kp1'],
            ptype,
            pot['a1'].value(),
            pot['a2'].value()
        )]

        if pot['rc'].value() > 0:
            lines[0] += " rc={:.2f}".format(pot['rc'].value())

        params = []
        if pot['uv'].value() != 0:
            params.append("            uv={:.2f} av={:.3f}".format(
                pot['uv'].value(), pot['av'].value()
            ))
            if pot['rv'].value() != 0:
                params[-1] += " rv={:.2f}".format(pot['rv'].value())

        if pot['uw'].value() != 0:
            params.append("            uw={:.2f}  aw={:.3f}".format(
                pot['uw'].value(), pot['aw'].value()
            ))
            if pot['rw'].value() != 0:
                params[-1] += " rw={:.2f}".format(pot['rw'].value())

        if pot['wd'].value() != 0:
            params.append("            wd={:.2f}  awd={:.3f}".format(
                pot['wd'].value(), pot['awd'].value()
            ))
            if pot['rwd'].value() != 0:
                params[-1] += " rwd={:.2f}".format(pot['rwd'].value())

        lines.extend(params)
        lines.append("           /\n\n")

        return '\n'.join(lines)

    def load_from_file(self, file_path):
        """Load input from a file"""
        with open(file_path, 'r') as f:
            content = f.read()

        # Parse namelists using regex
        self.parse_global(content)
        self.parse_system(content)
        self.parse_outgoing(content)
        self.parse_potentials(content)

    def parse_global(self, content):
        """Parse GLOBAL namelist"""
        global_match = re.search(r'&GLOBAL\s+(.*?)\s*/', content, re.DOTALL)
        if global_match:
            params = global_match.group(1)
            self.set_value(self.hcm, params, 'hcm', float)
            self.set_value(self.rmax, params, 'rmax', float)
            self.set_value(self.lmin, params, 'lmin', int)
            self.set_value(self.lmax, params, 'lmax', int)
            self.set_value(self.lxmax, params, 'lxmax', int)
            self.set_value(self.elab, params, 'elab', float)
            self.set_value(self.thmin, params, 'thmin', float)
            self.set_value(self.thmax, params, 'thmax', float)
            self.set_value(self.thinc, params, 'thinc', float)
            self.set_value(self.nx, params, 'nx', int)
            self.set_value(self.nr, params, 'nr', int)

            dwba_match = re.search(r'dwba=(\d+)', params)
            if dwba_match:
                self.dwba.setCurrentIndex(int(dwba_match.group(1)) - 1)

            printf_match = re.search(r'printf=([tf])', params)
            if printf_match:
                self.printf.setCurrentIndex(1 if printf_match.group(1) == 't' else 0)

    def parse_system(self, content):
        """Parse SYSTEM namelist"""
        system_match = re.search(r'&SYSTEM\s+(.*?)\s*/', content, re.DOTALL)
        if system_match:
            params = system_match.group(1)

            self.set_text(self.namep, params, 'namep')
            self.set_value(self.massp, params, 'massp', float)
            self.set_value(self.zp, params, 'zp', float)
            self.set_value(self.jp, params, 'jp', float)

            self.set_text(self.namet, params, 'namet')
            self.set_value(self.masst, params, 'masst', float)
            self.set_value(self.zt, params, 'zt', float)
            self.set_value(self.jt, params, 'jt', float)

            self.set_text(self.nameb, params, 'nameb')
            self.set_value(self.massb, params, 'massb', float)
            self.set_value(self.zb, params, 'zb', float)
            self.set_value(self.jb, params, 'jb', float)

            self.set_text(self.namex, params, 'namex')
            self.set_value(self.massx, params, 'massx', float)
            self.set_value(self.zx, params, 'zx', float)
            self.set_value(self.jx, params, 'jx', float)

            self.set_value(self.sbx, params, 'sbx', float)
            self.set_value(self.lbx, params, 'lbx', int)
            self.set_value(self.nodes, params, 'nodes', int)
            self.set_value(self.be, params, 'be', float)

    def parse_outgoing(self, content):
        """Parse OUTGOING namelist"""
        outgoing_match = re.search(r'&OUTGOING\s+(.*?)\s*/', content, re.DOTALL)
        if outgoing_match:
            params = outgoing_match.group(1)
            self.set_value(self.ecmbmin, params, 'ecmbmin', float)
            self.set_value(self.ecmbmax, params, 'ecmbmax', float)
            self.set_value(self.ecmbh, params, 'ecmbh', float)

    def parse_potentials(self, content):
        """Parse POTENTIAL namelists"""
        potential_matches = re.finditer(r'&POTENTIAL\s+(.*?)\s*/', content, re.DOTALL)

        pot_list = [self.pot1, self.pot2, self.pot3, self.pot4, self.pot5]
        for i, match in enumerate(potential_matches):
            if i >= len(pot_list):
                break

            params = match.group(1)
            if not params.strip():
                continue

            pot = pot_list[i]
            self.set_potential_values(pot, params)

    def set_potential_values(self, pot, params):
        """Set values for a potential from parsed parameters"""
        ptype_match = re.search(r'ptype=(\d+)', params)
        if ptype_match:
            ptype = int(ptype_match.group(1))
            if ptype <= 5:
                pot['ptype'].setCurrentIndex(ptype - 1)
            else:
                pot['ptype'].setCurrentIndex(5 + (ptype - 41))

        self.set_value(pot['a1'], params, 'a1', float)
        self.set_value(pot['a2'], params, 'a2', float)
        self.set_value(pot['rc'], params, 'rc', float)
        self.set_value(pot['uv'], params, 'uv', float)
        self.set_value(pot['av'], params, 'av', float)
        self.set_value(pot['rv'], params, 'rv', float)
        self.set_value(pot['uw'], params, 'uw', float)
        self.set_value(pot['aw'], params, 'aw', float)
        self.set_value(pot['rw'], params, 'rw', float)
        self.set_value(pot['wd'], params, 'wd', float)
        self.set_value(pot['awd'], params, 'awd', float)
        self.set_value(pot['rwd'], params, 'rwd', float)

    def set_value(self, widget, params, name, dtype):
        """Set widget value from parsed parameters"""
        match = re.search(rf'{name}=([\d.+-]+)', params)
        if match:
            try:
                value = dtype(match.group(1))
                widget.setValue(value)
            except:
                pass

    def set_text(self, widget, params, name):
        """Set text widget value from parsed parameters"""
        match = re.search(rf"{name}='([^']+)'", params)
        if match:
            widget.setText(match.group(1))

    def reset_all(self):
        """Reset all inputs to default values"""
        # Would need to recreate all widgets or set to defaults
        # For simplicity, just reset key values
        self.elab.setValue(25.5)
        self.namep.setText("d")
        self.namet.setText("93Nb")
        self.nameb.setText("p")
        self.namex.setText("n")
