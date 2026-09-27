"""Graphical interface (Qt / PySide6)."""

from __future__ import annotations

import os
import sys
import time
from dataclasses import asdict

from PySide6.QtCore import QObject, QSettings, Qt, QThread, QUrl, Signal
from PySide6.QtGui import QAction, QDesktopServices, QFont, QPixmap, QTextCursor
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QDoubleSpinBox, QFileDialog, QFormLayout, QGridLayout, QGroupBox,
    QHBoxLayout, QHeaderView, QLabel, QLineEdit, QListWidget, QMainWindow, QMessageBox,
    QPlainTextEdit, QProgressBar, QPushButton, QScrollArea, QSpinBox, QSplitter, QTableWidget,
    QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget,
)

from . import __version__
from .gui_import import ImportPanel
from .pipeline import Callbacks, Cancelled, Pipeline, Settings

STAGES = ["Import", "Kymographs", "Background", "Tracing", "Beat analysis"]


# ----------------------------------------------------------------------- worker

class Worker(QObject):
    log = Signal(str)
    progress = Signal(str, int, int)
    finished = Signal(object, str)   # report or None, error message

    def __init__(self, settings: Settings):
        super().__init__()
        self.settings = settings
        self._cancel = False

    def cancel(self):
        self._cancel = True

    def run(self):
        cb = Callbacks(log=self.log.emit, progress=self.progress.emit, cancelled=lambda: self._cancel)
        try:
            report = Pipeline(self.settings, cb).run()
            self.finished.emit(report, "")
        except Cancelled:
            self.finished.emit(None, "Cancelled")
        except Exception as e:  # show any problem to the user instead of crashing
            self.finished.emit(None, f"{type(e).__name__}: {e}")


# ----------------------------------------------------------------------- widgets

class PathRow(QWidget):
    def __init__(self, placeholder: str, pick_file: bool = False, file_filter: str = ""):
        super().__init__()
        self.pick_file, self.file_filter = pick_file, file_filter
        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        self.edit = QLineEdit()
        self.edit.setPlaceholderText(placeholder)
        btn = QPushButton("Browse…")
        btn.clicked.connect(self.browse)
        lay.addWidget(self.edit, 1)
        lay.addWidget(btn)

    def browse(self):
        start = self.edit.text() or os.path.expanduser("~")
        if self.pick_file:
            path, _ = QFileDialog.getOpenFileName(self, "Select file", start, self.file_filter)
        else:
            path = QFileDialog.getExistingDirectory(self, "Select folder", start)
        if path:
            self.edit.setText(path)

    def text(self) -> str:
        return self.edit.text().strip()

    def setText(self, t: str):
        self.edit.setText(t)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"tdtK Heart Analyzer {__version__}")
        self.resize(1100, 760)
        self.qs = QSettings("FlyHearts", "tdtK Analyzer")
        self.thread: QThread | None = None
        self.worker: Worker | None = None
        self.t_start = 0.0

        tabs = QTabWidget()
        tabs.addTab(self._run_tab(), "Analysis")
        self.import_panel = ImportPanel(self.settings, self._append_log)
        tabs.addTab(self.import_panel, "Import")
        tabs.addTab(self._settings_tab(), "Settings")
        tabs.addTab(self._results_tab(), "Results")
        self.tabs = tabs
        self.setCentralWidget(tabs)

        m = self.menuBar().addMenu("&Help")
        about = QAction("About", self)
        about.triggered.connect(self._about)
        m.addAction(about)
        self._load_settings()

    # ------------------------------------------------------------ tabs
    def _run_tab(self) -> QWidget:
        w = QWidget()
        v = QVBoxLayout(w)

        box = QGroupBox("Data")
        f = QFormLayout(box)
        self.movies = PathRow("Folder with the movies (.cxd, .nd2, .czi, .lif, .tif, …; sub-folders included)")
        self.output = PathRow("Output folder (created if needed)")
        self.mappings = PathRow("mappings.xlsx or mappings.csv", True, "Mappings (*.xlsx *.xls *.csv)")
        f.addRow("Movie folder:", self.movies)
        f.addRow("Output folder:", self.output)
        f.addRow("Genotype mappings:", self.mappings)
        v.addWidget(box)

        steps = QGroupBox("Steps")
        g = QGridLayout(steps)
        self.step1 = QCheckBox("1  Import movies, kymographs and beat direction  (check the movies on the Import tab)")
        self.step2 = QCheckBox("2  Background subtraction, edge tracing and quality control")
        self.step3 = QCheckBox("3  Beat analysis and summary tables")
        for i, cb in enumerate((self.step1, self.step2, self.step3)):
            cb.setChecked(True)
            g.addWidget(cb, i, 0)
        hint = QLabel("Files that already exist are not recomputed, so an interrupted run can simply be restarted.")
        hint.setStyleSheet("color: gray")
        g.addWidget(hint, 3, 0)
        v.addWidget(steps)

        prog = QGroupBox("Progress")
        pg = QGridLayout(prog)
        self.bars = {}
        for i, s in enumerate(STAGES):
            bar = QProgressBar()
            bar.setFormat("%v / %m")
            bar.setValue(0)
            bar.setMaximum(1)
            self.bars[s] = bar
            pg.addWidget(QLabel(s), i, 0)
            pg.addWidget(bar, i, 1)
        self.status = QLabel("Ready.")
        pg.addWidget(self.status, len(STAGES), 0, 1, 2)
        v.addWidget(prog)

        row = QHBoxLayout()
        self.run_btn = QPushButton("Run analysis")
        self.run_btn.setDefault(True)
        self.run_btn.clicked.connect(self.start)
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setEnabled(False)
        self.cancel_btn.clicked.connect(self.cancel)
        open_btn = QPushButton("Open output folder")
        open_btn.clicked.connect(lambda: self._open(self.output.text()))
        row.addWidget(self.run_btn)
        row.addWidget(self.cancel_btn)
        row.addStretch(1)
        row.addWidget(open_btn)
        v.addLayout(row)

        self.log_view = QPlainTextEdit()
        self.log_view.setReadOnly(True)
        self.log_view.setFont(QFont("Monospace", 9))
        self.log_view.setMaximumBlockCount(20000)
        v.addWidget(self.log_view, 1)
        return w

    def _settings_tab(self) -> QWidget:
        w = QWidget()
        f = QFormLayout(w)
        cpu = os.cpu_count() or 2
        self.workers = QSpinBox()
        self.workers.setRange(1, max(1, cpu))
        self.movie_workers = QSpinBox()
        self.movie_workers.setRange(1, max(1, cpu))
        self.min_size = QDoubleSpinBox()
        self.min_size.setRange(0, 100000)
        self.min_size.setSuffix(" MB")
        self.min_frames = QSpinBox()
        self.min_frames.setRange(1, 10_000_000)
        self.max_interval = QDoubleSpinBox()
        self.max_interval.setRange(0.1, 10000)
        self.max_interval.setDecimals(2)
        self.max_interval.setSuffix(" ms")
        self.ball = QDoubleSpinBox()
        self.ball.setRange(1, 1000)
        self.ball.setSuffix(" px")
        f.addRow("Parallel workers (tracing, analysis):", self.workers)
        f.addRow("Movies processed at the same time:", self.movie_workers)
        f.addRow(QLabel("Each movie is loaded completely into memory - increase only with enough RAM."))
        f.addRow("Skip .cxd movies smaller than:", self.min_size)
        f.addRow("Skip movies with fewer frames than:", self.min_frames)
        f.addRow("Skip movies with a frame interval above:", self.max_interval)
        f.addRow("Background rolling-ball radius:", self.ball)
        reset = QPushButton("Restore defaults")
        reset.clicked.connect(lambda: self._apply(Settings(movie_dir=self.movies.text(), output_dir=self.output.text(),
                                                            mappings_file=self.mappings.text())))
        f.addRow(reset)
        return w

    def _results_tab(self) -> QWidget:
        w = QWidget()
        v = QVBoxLayout(w)
        self.timing = QTableWidget(0, 3)
        self.timing.setHorizontalHeaderLabels(["Step", "Items", "Seconds"])
        self.timing.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.timing.setMaximumHeight(170)
        v.addWidget(QLabel("Timing of the last run (also saved as timing.csv):"))
        v.addWidget(self.timing)

        split = QSplitter(Qt.Horizontal)
        left = QWidget()
        lv = QVBoxLayout(left)
        lv.setContentsMargins(0, 0, 0, 0)
        lv.addWidget(QLabel("Summary tables and traced kymographs:"))
        self.files = QListWidget()
        self.files.itemDoubleClicked.connect(lambda it: self._open(it.data(Qt.UserRole)))
        self.files.currentItemChanged.connect(self._preview)
        lv.addWidget(self.files)
        refresh = QPushButton("Refresh")
        refresh.clicked.connect(self.refresh_results)
        lv.addWidget(refresh)
        split.addWidget(left)
        self.image = QLabel("Select a traced kymograph (*_traced.jpg) to preview it.\nDouble-click opens a file.")
        self.image.setAlignment(Qt.AlignCenter)
        scroll = QScrollArea()
        scroll.setWidget(self.image)
        scroll.setWidgetResizable(True)
        split.addWidget(scroll)
        split.setSizes([350, 700])
        v.addWidget(split, 1)
        return w

    # ------------------------------------------------------------ settings
    def settings(self) -> Settings:
        return Settings(movie_dir=self.movies.text(), output_dir=self.output.text(),
                        mappings_file=self.mappings.text(), run_kymographs=self.step1.isChecked(),
                        run_tracing=self.step2.isChecked(), run_analysis=self.step3.isChecked(),
                        workers=self.workers.value(), movie_workers=self.movie_workers.value(),
                        min_file_size_mb=self.min_size.value(), min_frames=self.min_frames.value(),
                        max_frame_interval_ms=self.max_interval.value(),
                        rolling_ball_radius=self.ball.value(), **self._import_settings())

    def _import_settings(self) -> dict:
        if not hasattr(self, "import_panel"):
            d = Settings()
            return {"channel": d.channel, "z_plane": d.z_plane, "rotate": d.rotate,
                    "default_interval_ms": d.default_interval_ms, "default_pixel_um": d.default_pixel_um}
        o = self.import_panel.options()
        return {"channel": o.channel, "z_plane": o.z_plane, "rotate": o.rotate,
                "default_interval_ms": o.default_interval_ms, "default_pixel_um": o.default_pixel_um}

    def _apply(self, s: Settings):
        self.movies.setText(s.movie_dir)
        self.output.setText(s.output_dir)
        self.mappings.setText(s.mappings_file)
        self.step1.setChecked(s.run_kymographs)
        self.step2.setChecked(s.run_tracing)
        self.step3.setChecked(s.run_analysis)
        self.workers.setValue(s.workers)
        self.movie_workers.setValue(s.movie_workers)
        self.min_size.setValue(s.min_file_size_mb)
        self.min_frames.setValue(s.min_frames)
        self.max_interval.setValue(s.max_frame_interval_ms)
        self.ball.setValue(s.rolling_ball_radius)
        if hasattr(self, "import_panel"):
            self.import_panel.set_options(s)

    def _load_settings(self):
        d = asdict(Settings())
        for k, default in d.items():
            val = self.qs.value(k, default)
            if isinstance(default, bool):
                val = val in (True, "true", "True", 1, "1")
            elif isinstance(default, int):
                val = int(val)
            elif isinstance(default, float):
                val = float(val)
            d[k] = val
        self._apply(Settings(**d))

    def _save_settings(self):
        for k, v in asdict(self.settings()).items():
            self.qs.setValue(k, v)

    # ------------------------------------------------------------ run
    def start(self):
        s = self.settings()
        problems = Pipeline(s).validate()
        if problems:
            QMessageBox.warning(self, "Cannot start", "\n".join(problems))
            return
        self._save_settings()
        if self.import_panel.infos and s.run_kymographs:
            self.import_panel.save()          # the run uses the choices made on the Import tab
        for bar in self.bars.values():
            bar.setMaximum(1)
            bar.setValue(0)
        self.log_view.clear()
        self.run_btn.setEnabled(False)
        self.cancel_btn.setEnabled(True)
        self.status.setText("Running…")
        self.t_start = time.time()
        self.thread = QThread()
        self.worker = Worker(s)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.log.connect(self._append_log)
        self.worker.progress.connect(self._on_progress)
        self.worker.finished.connect(self._on_finished)
        self.worker.finished.connect(self.thread.quit)
        self.thread.start()

    def cancel(self):
        if self.worker:
            self.worker.cancel()
            self.status.setText("Cancelling after the files that are being processed…")
            self.cancel_btn.setEnabled(False)

    def _append_log(self, line: str):
        self.log_view.appendPlainText(line)
        self.log_view.moveCursor(QTextCursor.End)

    def _on_progress(self, stage: str, done: int, total: int):
        bar = self.bars.get(stage)
        if bar:
            bar.setMaximum(max(total, 1))
            bar.setValue(done)
        elapsed = time.time() - self.t_start
        self.status.setText(f"{stage}: {done} of {total}   (elapsed {elapsed:.0f} s)")

    def _on_finished(self, report, error: str):
        self.run_btn.setEnabled(True)
        self.cancel_btn.setEnabled(False)
        if report is None:
            self.status.setText(error or "Stopped.")
            if error and error != "Cancelled":
                QMessageBox.critical(self, "Analysis stopped", error)
            return
        self.status.setText(f"Finished in {time.time() - self.t_start:.0f} s."
                            + (f"  {len(report.errors)} file(s) had problems - see the log." if report.errors else ""))
        self.timing.setRowCount(0)
        for t in report.timing:
            r = self.timing.rowCount()
            self.timing.insertRow(r)
            self.timing.setItem(r, 0, QTableWidgetItem(t.step))
            self.timing.setItem(r, 1, QTableWidgetItem(str(t.items) if t.items else ""))
            self.timing.setItem(r, 2, QTableWidgetItem(f"{t.seconds:.1f}"))
        self.refresh_results()
        self.tabs.setCurrentIndex(3)

    # ------------------------------------------------------------ results
    def refresh_results(self):
        self.files.clear()
        ex = os.path.join(self.output.text(), "balled", "excellent traces")
        if not os.path.isdir(ex):
            return
        names = sorted(os.listdir(ex))
        for n in [n for n in names if n.endswith(".csv") and not n.endswith(".tiff.csv")] + \
                 [n for n in names if n.endswith("_traced.jpg")]:
            self.files.addItem(n)
            it = self.files.item(self.files.count() - 1)
            it.setData(Qt.UserRole, os.path.join(ex, n))

    def _preview(self, item, _prev=None):
        if not item:
            return
        path = item.data(Qt.UserRole)
        if path.endswith(".jpg"):
            pm = QPixmap(path)
            if not pm.isNull():
                self.image.setPixmap(pm.scaled(max(pm.width(), 900), max(pm.height() * 4, 200),
                                               Qt.KeepAspectRatio, Qt.SmoothTransformation))
                return
        self.image.setText(os.path.basename(path) + "\n\nDouble-click to open.")

    def _open(self, path: str):
        if path and os.path.exists(path):
            QDesktopServices.openUrl(QUrl.fromLocalFile(path))

    def _about(self):
        QMessageBox.about(self, "tdtK Heart Analyzer",
                          f"<b>tdtK Heart Analyzer {__version__}</b><br><br>"
                          "Standalone analysis of tdTomato fly heart movies (.cxd).<br>"
                          "Python port of tdtK_Full_Analysis_script (R) - no R, Fiji or Java required.")

    def closeEvent(self, ev):
        if self.thread and self.thread.isRunning():
            if QMessageBox.question(self, "Quit", "An analysis is running. Cancel it and quit?") != QMessageBox.Yes:
                ev.ignore()
                return
            self.worker.cancel()
            self.thread.quit()
            self.thread.wait(5000)
        self._save_settings()
        ev.accept()


def main() -> int:
    app = QApplication.instance() or QApplication(sys.argv)
    app.setApplicationName("tdtK Heart Analyzer")
    w = MainWindow()
    w.show()
    return app.exec()
