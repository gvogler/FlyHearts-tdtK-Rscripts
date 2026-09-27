"""'Import' tab: find movies of any supported format, check and edit their metadata, preview them."""

from __future__ import annotations

import os
from typing import Callable

import numpy as np
from PySide6.QtCore import QObject, Qt, QThread, Signal
from PySide6.QtGui import QBrush, QColor, QGuiApplication, QImage, QPixmap
from PySide6.QtWidgets import (
    QAbstractItemView, QComboBox, QDoubleSpinBox, QFormLayout, QGroupBox, QHBoxLayout, QHeaderView, QLabel,
    QMessageBox, QProgressBar, QPushButton, QSplitter, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from .movies import (FORMAT_NAMES, ImportOptions, MovieInfo, apply_table, load_movie, load_table, movie_status,
                     save_table, scan_movies)
from .movies.bioformats_reader import available as bioformats_available

COLS = ["Use", "File", "Series", "Output name", "Format", "Frames", "Size (H×W)", "Channels",
        "Channel", "Rotate", "Frame interval (ms)", "Pixel size (µm)", "Status"]
C_USE, C_FILE, C_SERIES, C_NAME, C_FMT, C_T, C_SIZE, C_CH, C_CHSEL, C_ROT, C_INT, C_PIX, C_STATUS = range(len(COLS))
EDITABLE = {C_NAME, C_CHSEL, C_ROT, C_INT, C_PIX}
COLORS = {"ok": "#d9f2d9", "warn": "#fff3c4", "skip": "#e6e6e6", "error": "#f8d0d0"}


class ScanWorker(QObject):
    progress = Signal(int, int, str)
    finished = Signal(object, str)

    def __init__(self, folder: str, output: str):
        super().__init__()
        self.folder, self.output = folder, output

    def run(self):
        try:
            infos = scan_movies(self.folder, self.output or None, self.progress.emit)
            table = load_table(self.output) if self.output else None
            if table is not None:
                apply_table(infos, table)
            self.finished.emit(infos, "")
        except Exception as e:
            self.finished.emit(None, f"{type(e).__name__}: {e}")


def _to_qimage(img: np.ndarray) -> QImage:
    a = img.astype(np.float32)
    lo, hi = np.percentile(a, 0.5), np.percentile(a, 99.5)
    a = np.clip((a - lo) / max(hi - lo, 1e-9), 0, 1)
    a8 = np.ascontiguousarray((a * 255).astype(np.uint8))
    q = QImage(a8.data, a8.shape[1], a8.shape[0], a8.strides[0], QImage.Format_Grayscale8)
    return q.copy()


class ImportPanel(QWidget):
    def __init__(self, get_settings: Callable, log: Callable[[str], None]):
        super().__init__()
        self.get_settings = get_settings
        self.log = log
        self.infos: list[MovieInfo] = []
        self._updating = False
        v = QVBoxLayout(self)

        top = QHBoxLayout()
        self.scan_btn = QPushButton("Scan movie folder")
        self.scan_btn.clicked.connect(self.scan)
        self.save_btn = QPushButton("Save import table")
        self.save_btn.clicked.connect(self.save)
        self.bar = QProgressBar()
        self.bar.setMaximum(1)
        self.bar.setFormat("%v / %m files")
        self.info_lbl = QLabel("Scan the movie folder to see which movies will be analyzed.")
        top.addWidget(self.scan_btn)
        top.addWidget(self.save_btn)
        top.addWidget(self.bar, 1)
        v.addLayout(top)
        v.addWidget(self.info_lbl)

        defaults = QGroupBox("Defaults for all movies (a value in the table overrides them)")
        f = QHBoxLayout(defaults)
        form1, form2 = QFormLayout(), QFormLayout()
        self.channel = QComboBox()
        self.channel.setEditable(True)
        self.channel.addItems(["auto", "0", "1", "2", "3"])
        self.channel.setToolTip("'auto' picks a channel named like tdTomato/mCherry/RFP/561…, otherwise the first.\n"
                                "Enter a 0-based index or part of a channel name.")
        self.zplane = QComboBox()
        self.zplane.setEditable(True)
        self.zplane.addItems(["0", "max"])
        self.zplane.setToolTip("Z plane to use for recordings with several planes (index, or 'max' projection)")
        self.rotate = QComboBox()
        self.rotate.addItems(["0", "90", "180", "270", "auto"])
        self.rotate.setToolTip("The analysis expects the heart to run left-right.\n"
                               "Rotate counter-clockwise by 90/180/270°, or 'auto' to detect it.")
        self.def_int = QDoubleSpinBox()
        self.def_int.setRange(0, 10000)
        self.def_int.setDecimals(3)
        self.def_int.setSpecialValueText("from file")
        self.def_int.setSuffix(" ms")
        self.def_pix = QDoubleSpinBox()
        self.def_pix.setRange(0, 1000)
        self.def_pix.setDecimals(4)
        self.def_pix.setSpecialValueText("from file")
        self.def_pix.setSuffix(" µm")
        form1.addRow("Channel:", self.channel)
        form1.addRow("Z plane:", self.zplane)
        form1.addRow("Rotate:", self.rotate)
        form2.addRow("Frame interval if missing:", self.def_int)
        form2.addRow("Pixel size if missing:", self.def_pix)
        for w in (self.channel, self.zplane, self.rotate):
            w.currentTextChanged.connect(lambda *_: self.refresh_status())
        for w in (self.def_int, self.def_pix):
            w.valueChanged.connect(lambda *_: self.refresh_status())
        f.addLayout(form1)
        f.addLayout(form2)
        fmts = ", ".join(v for k, v in FORMAT_NAMES.items() if k != "bioformats")
        extra = ("Bio-Formats found: all other Bio-Formats formats are read too."
                 if bioformats_available() else "Other formats: install Bio-Formats 'bftools' (optional).")
        hint = QLabel(f"<small>Readable: {fmts}.<br>{extra}</small>")
        hint.setWordWrap(True)
        f.addWidget(hint, 1)
        v.addWidget(defaults)

        split = QSplitter(Qt.Vertical)
        self.table = QTableWidget(0, len(COLS))
        self.table.setHorizontalHeaderLabels(COLS)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(C_STATUS, QHeaderView.Stretch)
        self.table.itemChanged.connect(self._edited)
        self.table.currentCellChanged.connect(lambda r, *_: self.preview_btn.setEnabled(r >= 0))
        split.addWidget(self.table)

        pv = QWidget()
        pl = QVBoxLayout(pv)
        pl.setContentsMargins(0, 0, 0, 0)
        row = QHBoxLayout()
        self.preview_btn = QPushButton("Preview selected movie")
        self.preview_btn.setEnabled(False)
        self.preview_btn.clicked.connect(self.preview)
        row.addWidget(self.preview_btn)
        row.addWidget(QLabel("<small>Left: first frame. Right: movement (SD over time) - the beating heart "
                             "should appear as a bright band running <b>left-right</b>.</small>"), 1)
        pl.addLayout(row)
        imgs = QHBoxLayout()
        self.frame_lbl, self.sd_lbl = QLabel(), QLabel()
        for lbl in (self.frame_lbl, self.sd_lbl):
            lbl.setAlignment(Qt.AlignCenter)
            lbl.setMinimumHeight(140)
            lbl.setStyleSheet("background: #222; color: #ccc")
            imgs.addWidget(lbl, 1)
        pl.addLayout(imgs)
        split.addWidget(pv)
        split.setSizes([420, 220])
        v.addWidget(split, 1)
        self.thread = None

    # ------------------------------------------------------------ options
    def options(self) -> ImportOptions:
        return ImportOptions(channel=self.channel.currentText().strip() or "auto",
                             z_plane=self.zplane.currentText().strip() or "0", rotate=self.rotate.currentText(),
                             default_interval_ms=self.def_int.value(), default_pixel_um=self.def_pix.value())

    def set_options(self, s):
        self.channel.setCurrentText(s.channel)
        self.zplane.setCurrentText(s.z_plane)
        self.rotate.setCurrentText(s.rotate)
        self.def_int.setValue(s.default_interval_ms)
        self.def_pix.setValue(s.default_pixel_um)

    # ------------------------------------------------------------ scanning
    def scan(self):
        s = self.get_settings()
        if not os.path.isdir(s.movie_dir):
            QMessageBox.warning(self, "Scan", "Select the movie folder on the Analysis tab first.")
            return
        self.scan_btn.setEnabled(False)
        self.info_lbl.setText("Scanning…")
        self.thread = QThread()
        self.worker = ScanWorker(s.movie_dir, s.output_dir)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.progress.connect(self._scan_progress)
        self.worker.finished.connect(self._scanned)
        self.worker.finished.connect(self.thread.quit)
        self.thread.start()

    def _scan_progress(self, done, total, name):
        self.bar.setMaximum(max(total, 1))
        self.bar.setValue(done)
        if name:
            self.info_lbl.setText(f"Reading {name}")

    def _scanned(self, infos, error):
        self.scan_btn.setEnabled(True)
        if infos is None:
            self.info_lbl.setText("Scan failed.")
            QMessageBox.critical(self, "Scan failed", error)
            return
        self.infos = infos
        self.fill()
        self.log(f"Import: {len(infos)} movies found")

    def fill(self):
        self._updating = True
        self.table.setRowCount(len(self.infos))
        for r, m in enumerate(self.infos):
            use = QTableWidgetItem()
            use.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled | Qt.ItemIsSelectable)
            use.setCheckState(Qt.Checked if m.use else Qt.Unchecked)
            self.table.setItem(r, C_USE, use)
            chans = ", ".join(map(str, m.channel_names)) if m.channel_names else str(m.size_c)
            vals = {C_FILE: m.rel, C_SERIES: (m.series_name or str(m.series + 1)) if m.n_series > 1 else "",
                    C_NAME: m.name, C_FMT: m.format, C_T: str(m.size_t), C_SIZE: f"{m.size_y}×{m.size_x}",
                    C_CH: chans + (f"  (z: {m.size_z})" if m.size_z > 1 else ""),
                    C_CHSEL: m.channel or "", C_ROT: m.rotate or "",
                    C_INT: f"{m.time_interval * 1000:.4g}" if m.time_interval else "",
                    C_PIX: f"{m.pixel_size:.5g}" if m.pixel_size else ""}
            for c, text in vals.items():
                it = QTableWidgetItem(text)
                flags = Qt.ItemIsEnabled | Qt.ItemIsSelectable
                if c in EDITABLE:
                    flags |= Qt.ItemIsEditable
                    it.setToolTip("Double-click to edit" + (" (empty = default)" if c in (C_CHSEL, C_ROT) else ""))
                it.setFlags(flags)
                self.table.setItem(r, c, it)
            self.table.setItem(r, C_STATUS, QTableWidgetItem(""))
        self._updating = False
        self.refresh_status()

    def _edited(self, item):
        if self._updating:
            return
        r, c = item.row(), item.column()
        m = self.infos[r]
        text = item.text().strip()
        if c == C_USE:
            m.use = item.checkState() == Qt.Checked
        elif c == C_NAME and text:
            m.name = text
        elif c == C_CHSEL:
            m.channel = text or None
        elif c == C_ROT:
            m.rotate = text if text in ("0", "90", "180", "270", "auto") else None
        elif c in (C_INT, C_PIX):
            try:
                val = float(text.replace(",", ".")) if text else None
            except ValueError:
                val = None
            if c == C_INT:
                m.time_interval = val / 1000.0 if val else None
            else:
                m.pixel_size = val if val else None
        self.refresh_status()

    def refresh_status(self):
        s = self.get_settings()
        opts = self.options()
        counts = {"ok": 0, "warn": 0, "skip": 0, "error": 0}
        self._updating = True
        for r, m in enumerate(self.infos):
            level, msg = movie_status(m, opts, s.min_frames, s.max_frame_interval_ms / 1000.0)
            counts[level] += 1
            it = self.table.item(r, C_STATUS)
            if it is None:
                continue
            it.setText({"ok": "✓ ", "warn": "⚠ ", "skip": "– ", "error": "✗ "}[level] + msg)
            it.setToolTip(msg)
            for c in range(len(COLS)):
                cell = self.table.item(r, c)
                if cell:
                    cell.setBackground(QBrush(QColor(COLORS[level])))
        self._updating = False
        if self.infos:
            self.info_lbl.setText(f"{len(self.infos)} movies: {counts['ok'] + counts['warn']} ready, "
                                  f"{counts['skip']} skipped, {counts['error']} need attention (red).")

    def save(self) -> str | None:
        out = self.get_settings().output_dir
        if not self.infos:
            return None
        if not out:
            QMessageBox.warning(self, "Save", "Select an output folder on the Analysis tab first.")
            return None
        path = save_table(self.infos, out)
        self.log(f"Import table saved: {path}")
        return path

    # ------------------------------------------------------------ preview
    def preview(self):
        r = self.table.currentRow()
        if r < 0 or r >= len(self.infos):
            return
        m = self.infos[r]
        QGuiApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            movie, rot = load_movie(m, self.options(), max_frames=200)
        except Exception as e:
            QGuiApplication.restoreOverrideCursor()
            QMessageBox.warning(self, "Preview", f"{type(e).__name__}: {e}")
            return
        QGuiApplication.restoreOverrideCursor()
        sd = movie.astype(np.float32).std(axis=0)
        for lbl, img in ((self.frame_lbl, movie[0]), (self.sd_lbl, sd)):
            pm = QPixmap.fromImage(_to_qimage(img))
            lbl.setPixmap(pm.scaled(lbl.width() - 4, lbl.height() - 4, Qt.KeepAspectRatio, Qt.SmoothTransformation))
        self.info_lbl.setText(f"Preview of {m.name}: {movie.shape[0]} frames, rotated {rot}°, "
                              f"channel {m.channel_index(self.options())}")
