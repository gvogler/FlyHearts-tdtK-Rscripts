"""'Import' tab: find movies of any supported format, check their metadata, correct it, preview them.

* The table shows every movie; cells with missing (red) or doubtful (amber) metadata are
  highlighted, values entered by the user are shown in bold blue.
* The editor below the table shows where every value comes from and all problems found,
  and lets the user correct them for one movie or for all selected movies at once.
"""

from __future__ import annotations

import os
from typing import Callable

import numpy as np
from PySide6.QtCore import QObject, Qt, QThread, Signal
from PySide6.QtGui import QBrush, QColor, QFont, QGuiApplication, QImage, QPixmap
from PySide6.QtWidgets import (
    QAbstractItemView, QCheckBox, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFormLayout, QGridLayout,
    QGroupBox, QHBoxLayout, QHeaderView, QLabel, QLineEdit, QListWidget, QListWidgetItem, QMessageBox,
    QProgressBar, QPushButton, QSpinBox, QSplitter, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from .movies import (FORMAT_NAMES, ImportOptions, MovieInfo, apply_table, load_movie, load_table, movie_issues,
                     movie_status, save_table, scan_movies)
from .movies.bioformats_reader import available as bioformats_available
from .movies.metadata import format_date, parse_date, parse_interval, parse_pixel

COLS = ["Use", "File", "Series", "Output name", "Format", "Frames", "Size (H×W)", "Channels",
        "Channel", "Rotate", "Frame interval (ms)", "Pixel size (µm)", "Recorded", "Status"]
(C_USE, C_FILE, C_SERIES, C_NAME, C_FMT, C_T, C_SIZE, C_CH, C_CHSEL, C_ROT, C_INT, C_PIX, C_REC,
 C_STATUS) = range(len(COLS))
EDITABLE = {C_NAME, C_CHSEL, C_ROT, C_INT, C_PIX, C_REC}
FIELD_COL = {"time_interval": C_INT, "pixel_size": C_PIX, "created": C_REC, "channel": C_CHSEL, "axes": C_T,
             "file": C_FILE}
ROW_COLORS = {"ok": "#eaf7ea", "warn": "#fff8dc", "skip": "#eeeeee", "error": "#fbe3e3"}
CELL_COLORS = {"error": "#f4a6a6", "warn": "#ffd966"}
ICONS = {"error": "✗", "warn": "⚠", "info": "ℹ", "ok": "✓", "skip": "–"}


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
    return QImage(a8.data, a8.shape[1], a8.shape[0], a8.strides[0], QImage.Format_Grayscale8).copy()


def _fmt_interval(v):
    return f"{v * 1000:.5g}" if v else ""


def _fmt_pixel(v):
    return f"{v:.5g}" if v else ""


class PixelCalculator(QDialog):
    """pixel size = camera pixel × binning / (objective magnification × extra zoom)."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Pixel size calculator")
        f = QFormLayout(self)
        self.cam = QDoubleSpinBox(decimals=3, minimum=0.1, maximum=100, value=6.5, suffix=" µm")
        self.binning = QSpinBox(minimum=1, maximum=16, value=1)
        self.mag = QDoubleSpinBox(decimals=2, minimum=0.1, maximum=200, value=10, suffix=" ×")
        self.zoom = QDoubleSpinBox(decimals=3, minimum=0.01, maximum=10, value=1, suffix=" ×")
        self.result = QLabel()
        f.addRow("Camera pixel size:", self.cam)
        f.addRow("Binning:", self.binning)
        f.addRow("Objective magnification:", self.mag)
        f.addRow("Extra magnification (C-mount, optovar):", self.zoom)
        f.addRow("Pixel size in the image:", self.result)
        f.addRow(QLabel("<small>e.g. Hamamatsu ORCA-Flash4.0: 6.5 µm camera pixels</small>"))
        for w in (self.cam, self.binning, self.mag, self.zoom):
            w.valueChanged.connect(self._update)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        f.addRow(bb)
        self._update()

    def value(self) -> float:
        return self.cam.value() * self.binning.value() / (self.mag.value() * self.zoom.value())

    def _update(self):
        self.result.setText(f"<b>{self.value():.4g} µm</b>")


class MetadataEditor(QGroupBox):
    """Shows the metadata of the selected movie(s) with sources and problems; edits them."""
    applied = Signal()

    def __init__(self, panel: "ImportPanel"):
        super().__init__("Metadata of the selected movie")
        self.panel = panel
        g = QGridLayout(self)
        self.title = QLabel("Select a movie in the table.")
        self.title.setWordWrap(True)
        g.addWidget(self.title, 0, 0, 1, 4)

        self.interval = QLineEdit()
        self.interval.setPlaceholderText("e.g. 5 ms, 0.005 s or 200 fps")
        self.pixel = QLineEdit()
        self.pixel.setPlaceholderText("e.g. 0.65 or 650 nm")
        calc = QPushButton("Calculator…")
        calc.clicked.connect(self._calculator)
        self.recorded = QLineEdit()
        self.recorded.setPlaceholderText("YYYY-MM-DD HH:MM")
        self.channel = QComboBox()
        self.channel.setEditable(True)
        self.rotate = QComboBox()
        self.rotate.addItems(["(default)", "0", "90", "180", "270", "auto"])
        self.name = QLineEdit()
        self.use = QCheckBox("analyze this movie")
        self.src = {k: QLabel() for k in ("time_interval", "pixel_size", "created")}
        for lbl in self.src.values():
            lbl.setStyleSheet("color: gray")
        rows = [("Frame interval:", self.interval, self.src["time_interval"], None),
                ("Pixel size (µm):", self.pixel, self.src["pixel_size"], calc),
                ("Recorded:", self.recorded, self.src["created"], None),
                ("Channel:", self.channel, None, None), ("Rotate:", self.rotate, None, None),
                ("Output name:", self.name, None, None)]
        for i, (label, w, src, extra) in enumerate(rows, start=1):
            g.addWidget(QLabel(label), i, 0)
            g.addWidget(w, i, 1)
            if src is not None:
                g.addWidget(src, i, 2)
            if extra is not None:
                g.addWidget(extra, i, 3)
        g.addWidget(self.use, len(rows) + 1, 1)
        self.issues = QListWidget()
        self.issues.setMaximumHeight(110)
        g.addWidget(QLabel("Problems found:"), len(rows) + 2, 0, Qt.AlignTop)
        g.addWidget(self.issues, len(rows) + 2, 1, 1, 3)

        btns = QHBoxLayout()
        self.apply_btn = QPushButton("Apply")
        self.apply_btn.setToolTip("Save these values for the selected movie")
        self.apply_all_btn = QPushButton("Apply to all selected")
        self.apply_all_btn.setToolTip("Frame interval, pixel size, recording date, channel and rotation\n"
                                      "(only the fields you changed) for every selected movie")
        self.reset_btn = QPushButton("Reset to file values")
        self.apply_btn.clicked.connect(lambda: self._apply(all_selected=False))
        self.apply_all_btn.clicked.connect(lambda: self._apply(all_selected=True))
        self.reset_btn.clicked.connect(self._reset)
        for b in (self.apply_btn, self.apply_all_btn, self.reset_btn):
            btns.addWidget(b)
        btns.addStretch(1)
        g.addLayout(btns, len(rows) + 3, 0, 1, 4)
        g.setColumnStretch(1, 1)
        self.current: MovieInfo | None = None
        self._shown = {}
        self.setEnabled(False)

    # ------------------------------------------------------------ show
    def show_movie(self, m: MovieInfo | None, n_selected: int):
        self.current = m
        self.setEnabled(m is not None)
        if m is None:
            self.title.setText("Select a movie in the table.")
            self.issues.clear()
            return
        extra = f"   —   {n_selected} movies selected" if n_selected > 1 else ""
        chans = ", ".join(map(str, m.channel_names)) or str(m.size_c)
        self.title.setText(f"<b>{m.name}</b>  ({m.format}; {m.size_t} frames of {m.size_y}×{m.size_x} px; "
                           f"channels: {chans}; z planes: {m.size_z}){extra}")
        self.interval.setText(f"{m.time_interval * 1000:.6g} ms" if m.time_interval else "")
        self.pixel.setText(_fmt_pixel(m.pixel_size))
        self.recorded.setText(format_date(m.created_unix))
        for key, lbl in self.src.items():
            s = m.sources.get(key)
            orig = (m.original or {}).get({"created": "created_unix"}.get(key, key))
            if s == "entered" and orig:
                of = {"time_interval": lambda v: f"{v * 1000:.4g} ms", "pixel_size": lambda v: f"{v:.4g} µm",
                      "created": format_date}[key](orig)
                s = f"entered by you (file: {of})"
            elif s == "entered":
                s = "entered by you (not in the file)"
            lbl.setText(s or "missing")
            lbl.setStyleSheet("color: #b00000" if not s else ("color: #1a4fb3" if "entered" in s else "color: gray"))
        if m.time_interval:
            self.src["time_interval"].setText(self.src["time_interval"].text() + f"  ({1 / m.time_interval:.4g} fps)")
        self.channel.blockSignals(True)
        self.channel.clear()
        self.channel.addItem("(default)")
        for i in range(m.size_c):
            nm = m.channel_names[i] if i < len(m.channel_names) else ""
            self.channel.addItem(f"{i}" + (f": {nm}" if nm else ""))
        self.channel.setCurrentIndex(0)
        if m.channel not in (None, ""):
            self.channel.setCurrentText(str(m.channel))
        self.channel.blockSignals(False)
        self.rotate.setCurrentText(m.rotate if m.rotate else "(default)")
        self.name.setText(m.name)
        self.use.setChecked(m.use)
        self.issues.clear()
        issues = movie_issues(m, self.panel.options())
        for it in issues:
            li = QListWidgetItem(f"{ICONS[it.level]}  {it}")
            li.setForeground(QBrush(QColor({"error": "#b00000", "warn": "#8a6d00", "info": "#555555"}[it.level])))
            self.issues.addItem(li)
        if not issues:
            self.issues.addItem(QListWidgetItem("✓  metadata complete"))
        self._shown = self._values()

    def _values(self) -> dict:
        ch = self.channel.currentText().strip()
        ch = "" if ch == "(default)" else ch.split(":")[0].strip()
        rot = self.rotate.currentText()
        return {"interval": self.interval.text().strip(), "pixel": self.pixel.text().strip(),
                "recorded": self.recorded.text().strip(), "channel": ch,
                "rotate": "" if rot == "(default)" else rot, "name": self.name.text().strip(),
                "use": self.use.isChecked()}

    # ------------------------------------------------------------ edit
    def _calculator(self):
        d = PixelCalculator(self)
        if d.exec():
            self.pixel.setText(f"{d.value():.5g}")

    def _apply(self, all_selected: bool):
        if self.current is None:
            return
        v = self._values()
        try:
            interval = parse_interval(v["interval"]) if v["interval"] else None
            pixel = parse_pixel(v["pixel"]) if v["pixel"] else None
            recorded = parse_date(v["recorded"]) if v["recorded"] else None
        except ValueError as e:
            QMessageBox.warning(self, "Invalid value", str(e))
            return
        changed = {k for k in v if v[k] != self._shown.get(k)}
        targets = self.panel.selected_movies() if all_selected else [self.current]
        for m in targets:
            single = m is self.current
            if single or "interval" in changed:
                if interval != m.time_interval:
                    m.set_value("time_interval", interval)
            if single or "pixel" in changed:
                if pixel != m.pixel_size:
                    m.set_value("pixel_size", pixel)
            if single or "recorded" in changed:
                if recorded is None or format_date(recorded) != format_date(m.created_unix):
                    m.set_value("created_unix", recorded)
            if single or "channel" in changed:
                m.channel = v["channel"] or None
            if single or "rotate" in changed:
                m.rotate = v["rotate"] or None
            if single:
                m.name = v["name"] or m.name
                m.use = v["use"]
            elif "use" in changed:
                m.use = v["use"]
        self.applied.emit()

    def _reset(self):
        for m in self.panel.selected_movies():
            m.reset_to_file()
        self.applied.emit()


class ImportPanel(QWidget):
    def __init__(self, get_settings: Callable, log: Callable[[str], None]):
        super().__init__()
        self.get_settings = get_settings
        self.log = log
        self.infos: list[MovieInfo] = []
        self._updating = False
        self.thread = None
        v = QVBoxLayout(self)

        top = QHBoxLayout()
        self.scan_btn = QPushButton("Scan movie folder")
        self.scan_btn.clicked.connect(self.scan)
        self.save_btn = QPushButton("Save import table")
        self.save_btn.clicked.connect(self.save)
        self.only_problems = QCheckBox("show only movies that need attention")
        self.only_problems.toggled.connect(lambda *_: self.refresh_status())
        self.bar = QProgressBar()
        self.bar.setMaximum(1)
        self.bar.setFormat("%v / %m files")
        top.addWidget(self.scan_btn)
        top.addWidget(self.save_btn)
        top.addWidget(self.only_problems)
        top.addWidget(self.bar, 1)
        v.addLayout(top)
        self.info_lbl = QLabel("Scan the movie folder to see which movies will be analyzed.")
        v.addWidget(self.info_lbl)

        defaults = QGroupBox("Defaults for all movies (values of a movie override them)")
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
        self.rotate = QComboBox()
        self.rotate.addItems(["0", "90", "180", "270", "auto"])
        self.rotate.setToolTip("The analysis expects the heart to run left-right.")
        self.def_int = QDoubleSpinBox(decimals=3, maximum=10000, suffix=" ms")
        self.def_int.setSpecialValueText("none")
        self.def_pix = QDoubleSpinBox(decimals=4, maximum=1000, suffix=" µm")
        self.def_pix.setSpecialValueText("none")
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
        hint = QLabel(f"<small>Readable: {fmts}.<br>{extra}<br>Red cells: missing metadata. Amber: doubtful "
                      "values. <b style='color:#1a4fb3'>Blue</b>: entered by you. Double-click a cell to edit, "
                      "or use the editor below (select several rows to correct them together).</small>")
        hint.setWordWrap(True)
        f.addWidget(hint, 1)
        v.addWidget(defaults)

        split = QSplitter(Qt.Vertical)
        self.table = QTableWidget(0, len(COLS))
        self.table.setHorizontalHeaderLabels(COLS)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(C_STATUS, QHeaderView.Stretch)
        self.table.horizontalHeader().setMinimumSectionSize(40)
        self.table.setWordWrap(False)             # full status text in the tooltip and the editor below
        self.table.itemChanged.connect(self._edited)
        self.table.itemSelectionChanged.connect(self._selection)
        split.addWidget(self.table)

        bottom = QSplitter(Qt.Horizontal)
        self.editor = MetadataEditor(self)
        self.editor.applied.connect(self._after_edit)
        bottom.addWidget(self.editor)
        pv = QWidget()
        pl = QVBoxLayout(pv)
        pl.setContentsMargins(0, 0, 0, 0)
        row = QHBoxLayout()
        self.preview_btn = QPushButton("Preview")
        self.preview_btn.setEnabled(False)
        self.preview_btn.clicked.connect(self.preview)
        row.addWidget(self.preview_btn)
        row.addWidget(QLabel("<small>Top: first frame. Bottom: movement (SD over time) - the heart should be a "
                             "bright band running <b>left-right</b>.</small>"), 1)
        pl.addLayout(row)
        self.frame_lbl, self.sd_lbl = QLabel(), QLabel()
        for lbl in (self.frame_lbl, self.sd_lbl):
            lbl.setAlignment(Qt.AlignCenter)
            lbl.setMinimumHeight(110)
            lbl.setStyleSheet("background: #222; color: #ccc")
            pl.addWidget(lbl, 1)
        bottom.addWidget(pv)
        bottom.setSizes([620, 500])
        split.addWidget(bottom)
        split.setSizes([330, 380])
        v.addWidget(split, 1)

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

    def selected_movies(self) -> list[MovieInfo]:
        rows = sorted({i.row() for i in self.table.selectedIndexes()})
        return [self.infos[r] for r in rows if r < len(self.infos)]

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
            self.table.setItem(r, C_USE, use)
            for c in range(1, len(COLS)):
                it = QTableWidgetItem("")
                flags = Qt.ItemIsEnabled | Qt.ItemIsSelectable
                if c in EDITABLE:
                    flags |= Qt.ItemIsEditable
                it.setFlags(flags)
                self.table.setItem(r, c, it)
        self._updating = False
        self.refresh_status()
        if self.infos:
            self.table.selectRow(0)

    def _row_values(self, m: MovieInfo) -> dict:
        chans = ", ".join(map(str, m.channel_names)) if m.channel_names else str(m.size_c)
        return {C_FILE: m.rel, C_SERIES: (m.series_name or str(m.series + 1)) if m.n_series > 1 else "",
                C_NAME: m.name, C_FMT: m.format, C_T: str(m.size_t), C_SIZE: f"{m.size_y}×{m.size_x}",
                C_CH: chans + (f"  (z: {m.size_z})" if m.size_z > 1 else ""),
                C_CHSEL: m.channel or "", C_ROT: m.rotate or "", C_INT: _fmt_interval(m.time_interval),
                C_PIX: _fmt_pixel(m.pixel_size), C_REC: format_date(m.created_unix)}

    def refresh_status(self):
        s = self.get_settings()
        opts = self.options()
        counts = {"ok": 0, "warn": 0, "skip": 0, "error": 0}
        self._updating = True
        bold = QFont()
        bold.setBold(True)
        for r, m in enumerate(self.infos):
            level, msg = movie_status(m, opts, s.min_frames, s.max_frame_interval_ms / 1000.0)
            counts[level] += 1
            issues = movie_issues(m, opts)
            vals = self._row_values(m)
            self.table.item(r, C_USE).setCheckState(Qt.Checked if m.use else Qt.Unchecked)
            for c, text in vals.items():
                self.table.item(r, c).setText(text)
            st = self.table.item(r, C_STATUS)
            st.setText(f"{ICONS[level]} {msg}")
            st.setToolTip("\n".join(f"{ICONS[i.level]} {i}" for i in issues) or msg)
            for c in range(len(COLS)):
                cell = self.table.item(r, c)
                cell.setBackground(QBrush(QColor(ROW_COLORS[level])))
                cell.setForeground(QBrush(QColor("black")))
                cell.setFont(QFont())
                if c in (C_INT, C_PIX, C_REC):
                    key = {C_INT: "time_interval", C_PIX: "pixel_size", C_REC: "created"}[c]
                    cell.setToolTip(f"source: {m.sources.get(key, 'missing')}")
            for i in issues:                       # highlight the cells that have a problem
                col = FIELD_COL.get(i.field)
                if col is not None and i.level in CELL_COLORS:
                    cell = self.table.item(r, col)
                    cell.setBackground(QBrush(QColor(CELL_COLORS[i.level])))
                    cell.setToolTip(cell.toolTip() + f"\n{ICONS[i.level]} {i.message}")
            for key, col in (("time_interval", C_INT), ("pixel_size", C_PIX), ("created", C_REC)):
                if m.sources.get(key) == "entered":
                    cell = self.table.item(r, col)
                    cell.setForeground(QBrush(QColor("#1a4fb3")))
                    cell.setFont(bold)
            hide = self.only_problems.isChecked() and level in ("ok", "skip")
            self.table.setRowHidden(r, hide)
        self._updating = False
        if self.infos:
            self.info_lbl.setText(f"{len(self.infos)} movies: {counts['ok']} ready, {counts['warn']} ready but "
                                  f"please check (amber), {counts['error']} with missing metadata (red), "
                                  f"{counts['skip']} skipped.")
        self._selection()

    def _selection(self):
        sel = self.selected_movies()
        r = self.table.currentRow()
        cur = self.infos[r] if 0 <= r < len(self.infos) and self.infos[r] in sel else (sel[0] if sel else None)
        self.editor.show_movie(cur, len(sel))
        self.preview_btn.setEnabled(cur is not None)

    def _after_edit(self):
        self.refresh_status()
        self.save(quiet=True)

    def _edited(self, item):
        if self._updating:
            return
        r, c = item.row(), item.column()
        m = self.infos[r]
        text = item.text().strip()
        try:
            if c == C_USE:
                m.use = item.checkState() == Qt.Checked
            elif c == C_NAME and text:
                m.name = text
            elif c == C_CHSEL:
                m.channel = text or None
            elif c == C_ROT:
                m.rotate = text if text in ("0", "90", "180", "270", "auto") else None
            elif c == C_INT:
                v = parse_interval(text) if text else None
                if v != m.time_interval:
                    m.set_value("time_interval", v)
            elif c == C_PIX:
                v = parse_pixel(text) if text else None
                if v != m.pixel_size:
                    m.set_value("pixel_size", v)
            elif c == C_REC:
                v = parse_date(text) if text else None
                if format_date(v) != format_date(m.created_unix):
                    m.set_value("created_unix", v)
        except ValueError as e:
            QMessageBox.warning(self, "Invalid value", str(e))
        self._after_edit()

    def save(self, quiet: bool = False) -> str | None:
        out = self.get_settings().output_dir
        if not self.infos:
            return None
        if not out:
            if not quiet:
                QMessageBox.warning(self, "Save", "Select an output folder on the Analysis tab first.")
            return None
        path = save_table(self.infos, out)
        if not quiet:
            self.log(f"Import table saved: {path}")
        return path

    # ------------------------------------------------------------ preview
    def preview(self):
        m = self.editor.current
        if m is None:
            return
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
            lbl.setPixmap(pm.scaled(max(lbl.width() - 4, 50), max(lbl.height() - 4, 50), Qt.KeepAspectRatio,
                                    Qt.SmoothTransformation))
        self.info_lbl.setText(f"Preview of {m.name}: {movie.shape[0]} frames, rotated {rot}°, "
                              f"channel {m.channel_index(self.options())}")
