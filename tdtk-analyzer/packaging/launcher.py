"""Entry point for the PyInstaller build (GUI without arguments, CLI with arguments)."""

import sys

from tdtk_analyzer.__main__ import main

if __name__ == "__main__":
    sys.exit(main())
