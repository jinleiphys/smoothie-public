"""
Intelligent path detection utilities for SMOOTHIE GUI
"""

import os
import sys


def find_repo_root():
    """
    Find the SMOOTHIE repository root directory.

    Searches upward from the GUI script location to find the repo root.
    Returns the absolute path to the repository root.
    """
    # Start from the directory containing this script
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Navigate up to find repo root (should be one level up from smoothie_gui/)
    repo_root = os.path.dirname(current_dir)

    # Verify this looks like the SMOOTHIE repo by checking for key directories
    smoothie_dir = os.path.join(repo_root, "smoothie")
    cm2lab_dir = os.path.join(repo_root, "cm2lab")

    if os.path.isdir(smoothie_dir) and os.path.isdir(cm2lab_dir):
        return repo_root

    # Fallback: search upward more aggressively
    check_dir = current_dir
    for _ in range(5):  # Check up to 5 levels up
        check_dir = os.path.dirname(check_dir)
        smoothie_dir = os.path.join(check_dir, "smoothie")
        cm2lab_dir = os.path.join(check_dir, "cm2lab")

        if os.path.isdir(smoothie_dir) and os.path.isdir(cm2lab_dir):
            return check_dir

    # If we can't find it, return None and let caller handle it
    return None


def find_executable(name, repo_root=None):
    """
    Intelligently find an executable (smoothie or cm2lab).

    Search order:
    1. Environment variable (SMOOTHIE_EXE or CM2LAB_EXE)
    2. Relative to detected repo root
    3. In PATH
    4. Common installation locations

    Args:
        name: Name of executable ('smoothie' or 'cm2lab')
        repo_root: Repository root directory (auto-detected if None)

    Returns:
        Absolute path to executable if found, None otherwise
    """
    # Check environment variable first
    env_var = f"{name.upper()}_EXE"
    if env_var in os.environ:
        exe_path = os.environ[env_var]
        if os.path.isfile(exe_path) and os.access(exe_path, os.X_OK):
            return exe_path

    # Auto-detect repo root if not provided
    if repo_root is None:
        repo_root = find_repo_root()

    # Search in repo-relative location
    if repo_root:
        relative_path = os.path.join(repo_root, name, name)
        if os.path.isfile(relative_path) and os.access(relative_path, os.X_OK):
            return relative_path

    # Search in PATH
    for path_dir in os.environ.get("PATH", "").split(os.pathsep):
        exe_path = os.path.join(path_dir, name)
        if os.path.isfile(exe_path) and os.access(exe_path, os.X_OK):
            return exe_path

    # Common installation locations (Unix-like systems)
    common_locations = [
        f"/usr/local/bin/{name}",
        f"/opt/smoothie/bin/{name}",
        os.path.expanduser(f"~/bin/{name}"),
    ]

    for location in common_locations:
        if os.path.isfile(location) and os.access(location, os.X_OK):
            return location

    return None


def get_default_test_directory(repo_root=None):
    """
    Get the default test directory for SMOOTHIE input files.

    Args:
        repo_root: Repository root directory (auto-detected if None)

    Returns:
        Absolute path to test directory if found, user's home directory otherwise
    """
    if repo_root is None:
        repo_root = find_repo_root()

    if repo_root:
        test_dir = os.path.join(repo_root, "smoothie", "test")
        if os.path.isdir(test_dir):
            return test_dir

    # Fallback to user's home directory
    return os.path.expanduser("~")


def get_executable_info(name, repo_root=None):
    """
    Get information about an executable for user feedback.

    Returns a tuple of (path, status_message, is_found)
    """
    exe_path = find_executable(name, repo_root)

    if exe_path:
        return (exe_path, f"Found {name} at: {exe_path}", True)
    else:
        if repo_root:
            expected_path = os.path.join(repo_root, name, name)
            message = (
                f"{name.upper()} executable not found.\n\n"
                f"Expected location: {expected_path}\n\n"
                f"Please build {name.upper()} first:\n"
                f"  cd {os.path.join(repo_root, name)}\n"
                f"  make\n\n"
                f"Or set {name.upper()}_EXE environment variable to the executable path."
            )
        else:
            message = (
                f"{name.upper()} executable not found.\n\n"
                f"Could not detect repository location.\n"
                f"Please set {name.upper()}_EXE environment variable to the executable path."
            )
        return (None, message, False)


# Cache the repo root to avoid repeated filesystem searches
_cached_repo_root = None


def get_repo_root():
    """Get cached repo root or detect it"""
    global _cached_repo_root
    if _cached_repo_root is None:
        _cached_repo_root = find_repo_root()
    return _cached_repo_root
