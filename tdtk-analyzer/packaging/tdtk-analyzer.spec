# PyInstaller spec: builds a standalone "tdtK Heart Analyzer" app (no Python/R needed).
#   pyinstaller packaging/tdtk-analyzer.spec --noconfirm
# Output: dist/tdtK-Heart-Analyzer/ (Windows/Linux) or dist/tdtK Heart Analyzer.app (macOS)
import sys

from PyInstaller.utils.hooks import collect_all

block_cipher = None
datas, binaries, hidden = [], [], ["openpyxl", "olefile", "tifffile"]
# readers are imported dynamically (importlib) - list them explicitly
hidden += [f"tdtk_analyzer.movies.{k}_reader" for k in
           ("cxd", "tiff", "nd2", "czi", "lif", "oif", "video", "bioformats")]
for pkg in ("nd2", "readlif", "oiffile", "pylibCZIrw", "imagecodecs", "av"):
    d, b, h = collect_all(pkg)
    datas, binaries, hidden = datas + d, binaries + b, hidden + h
excludes = ["cryptography", "tkinter", "matplotlib", "IPython", "pytest", "PySide6.Qt3DCore", "PySide6.QtWebEngineCore",
            "PySide6.QtWebEngineWidgets", "PySide6.QtQuick", "PySide6.QtQml", "PySide6.QtMultimedia",
            "PySide6.QtCharts", "PySide6.QtDataVisualization", "PySide6.QtPdf"]

a = Analysis(["launcher.py"], pathex=[".."], hiddenimports=hidden, datas=datas, binaries=binaries,
             excludes=excludes, cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="tdtK-Heart-Analyzer",
          console=False, upx=False)
coll = COLLECT(exe, a.binaries, a.zipfiles, a.datas, name="tdtK-Heart-Analyzer", upx=False)
if sys.platform == "darwin":
    app = BUNDLE(coll, name="tdtK Heart Analyzer.app", bundle_identifier="org.flyhearts.tdtk-analyzer")
