# PyInstaller spec: builds a standalone "tdtK Heart Analyzer" app (no Python/R needed).
#   pyinstaller packaging/tdtk-analyzer.spec --noconfirm
# Output: dist/tdtK-Heart-Analyzer/ (Windows/Linux) or dist/tdtK Heart Analyzer.app (macOS)
import sys

block_cipher = None
excludes = ["tkinter", "matplotlib", "IPython", "pytest", "PySide6.Qt3DCore", "PySide6.QtWebEngineCore",
            "PySide6.QtWebEngineWidgets", "PySide6.QtQuick", "PySide6.QtQml", "PySide6.QtMultimedia",
            "PySide6.QtCharts", "PySide6.QtDataVisualization", "PySide6.QtPdf"]

a = Analysis(["launcher.py"], pathex=[".."], hiddenimports=["openpyxl", "olefile", "tifffile"],
             excludes=excludes, cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="tdtK-Heart-Analyzer",
          console=False, upx=False)
coll = COLLECT(exe, a.binaries, a.zipfiles, a.datas, name="tdtK-Heart-Analyzer", upx=False)
if sys.platform == "darwin":
    app = BUNDLE(coll, name="tdtK Heart Analyzer.app", bundle_identifier="org.flyhearts.tdtk-analyzer")
