"""CPU compatibility checks used before importing native application dependencies."""

import ctypes
import platform
import subprocess


def supports_avx2():
    """Return whether the current CPU exposes AVX2 instructions."""
    system = platform.system()

    if system == "Windows":
        return bool(ctypes.windll.kernel32.IsProcessorFeaturePresent(40))

    if system == "Linux":
        try:
            with open("/proc/cpuinfo", encoding="ascii") as cpuinfo:
                return any("avx2" in line.lower() for line in cpuinfo)
        except (OSError, UnicodeError):
            return True

    if system == "Darwin":
        try:
            result = subprocess.run(
                ["sysctl", "-n", "hw.optional.avx2_0"],
                capture_output=True,
                check=False,
                text=True,
            )
            return result.stdout.strip() == "1"
        except OSError:
            return True

    return True


def show_compatibility_warning():
    """Show the fallback download location without importing the app runtime."""
    message = (
        "This computer does not support the CPU instructions required by the "
        "standard GNNPCSAFT WEBAPP build.\n\n"
        "Please download the CPU-compatible build."
    )

    if platform.system() == "Windows":
        ctypes.windll.user32.MessageBoxW(
            None,
            message,
            "GNNPCSAFT WEBAPP - incompatible CPU",
            0x10,
        )
    else:
        print(message)
