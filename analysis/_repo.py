"""Repository paths, so the analysis scripts run from any working directory.

Every script here reads from data/ and structures/ at the repository root.
Importing this rather than hardcoding a path is what lets someone clone the
repo and run the scripts without editing them first.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
STRUCTURES = ROOT / "structures"
SURVIVORS = STRUCTURES / "survivors"
TARGETS = ROOT / "targets"

# Scripts that predate this module build filenames by string concatenation,
# e.g. open(B + "liu_metrics.tsv"). Keeping a trailing separator lets them
# carry on doing that unchanged.
B = str(DATA) + "/"


def sibling(name):
    """Absolute path to another script in this directory.

    Several scripts reuse a helper by exec'ing a sibling file. That only works
    when the interpreter's working directory happens to be analysis/, which it
    usually is not.
    """
    return str(Path(__file__).resolve().parent / name)
