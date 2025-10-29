"""
SMOOTHIE runner for executing calculations
"""

from PySide6.QtCore import QObject, QProcess, Signal
import os


class SmoothieRunner(QObject):
    """Class to run SMOOTHIE calculations"""

    output_ready = Signal(str)
    error_ready = Signal(str)
    finished = Signal(int)
    started = Signal()

    def __init__(self):
        super().__init__()
        self.process = None

    def run(self, executable, input_file, working_dir=None):
        """Run SMOOTHIE with the given input file in specified working directory"""
        if self.process is not None and self.process.state() == QProcess.Running:
            self.error_ready.emit("A calculation is already running")
            return

        self.process = QProcess()
        self.process.setProcessChannelMode(QProcess.MergedChannels)

        # Connect signals
        self.process.readyReadStandardOutput.connect(self.handle_output)
        self.process.readyReadStandardError.connect(self.handle_error)
        self.process.finished.connect(self.handle_finished)
        self.process.started.connect(self.handle_started)

        # Set working directory (use specified dir or default to input file's dir)
        if working_dir is None:
            working_dir = os.path.dirname(input_file)
        self.process.setWorkingDirectory(working_dir)

        # Start the process with input redirection
        self.process.start(executable, [])

        # Write input file content to stdin
        if os.path.exists(input_file):
            with open(input_file, 'r') as f:
                input_content = f.read()
            self.process.write(input_content.encode())
            self.process.closeWriteChannel()

    def stop(self):
        """Stop the running calculation"""
        if self.process is not None and self.process.state() == QProcess.Running:
            self.process.kill()

    def handle_output(self):
        """Handle standard output from the process"""
        if self.process:
            data = self.process.readAllStandardOutput()
            text = bytes(data).decode('utf-8', errors='replace')
            self.output_ready.emit(text.strip())

    def handle_error(self):
        """Handle standard error from the process"""
        if self.process:
            data = self.process.readAllStandardError()
            text = bytes(data).decode('utf-8', errors='replace')
            self.error_ready.emit(text.strip())

    def handle_finished(self, exit_code, exit_status):
        """Handle process completion"""
        self.finished.emit(exit_code)

    def handle_started(self):
        """Handle process start"""
        self.started.emit()
