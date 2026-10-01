"""'Review' tab: inspect traced kymographs (M-modes) and choose which ones go into the analysis.

Traces the automatic quality control rejected ('good traces', 'bad traces') can be
added to 'excellent traces', and accepted ones can be removed. Decisions are kept
in balled/manual_curation.csv (see curation.py) and survive a new QC run.
"""

from __future__ import annotations

import os
from typing import Callable

import numpy as np
from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QBrush, QColor, QFont, QKeySequence, QPainter, QPen, QPixmap, QPolygonF, QShortcut
from PySide6.QtWidgets import (
    QAbstractItemView, QComboBox, QHBoxLayout, QHeaderView, QLabel, QMessageBox, QPushButton, QScrollArea,
    QSplitter, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

from . import curation, rio

COLS = ["In analysis", "Kymograph", "Movie", "Xpos", "Automatic QC", "Your decision",
        "SD", "Autocorr. up", "Autocorr. down", "Dist. corr."]
C_IN, C_NAME, C_MOVIE, C_X, C_AUTO, C_DEC, C_SD, C_AU, C_AD, C_DC = range(len(COLS))
AUTO_TEXT = {"excellent": "excellent", "rescued": "good (rescued*)", "good": "good", "bad": "bad"}
SHOW = {"Good traces": lambda t: t["automatic"].isin(["good", "rescued"]),
        "Bad traces": lambda t: t["automatic"] == "bad",
        "In the analysis": lambda t: t["in_analysis"],
        "Not in the analysis": lambda t: ~t["in_analysis"],
        "Your decisions": lambda t: t["decision"] != "",
        "All traces": lambda t: t["automatic"].notna()}


def diameter_plot(csv_path: str, width: int, height: int) -> QPixmap:
    """Heart diameter over time from a trace CSV ('<...>.tiff.csv')."""
    pm = QPixmap(max(width, 200), max(height, 80))
    pm.fill(QColor("white"))
    try:
        d = rio.read_r_csv(csv_path)
        t = d.iloc[:, 0].to_numpy(float)
        y = d.iloc[:, 1].to_numpy(float)
    except Exception:
        return pm
    ok = np.isfinite(t) & np.isfinite(y)
    t, y = t[ok], y[ok]
    if t.size < 2:
        return pm
    p = QPainter(pm)
    p.setRenderHint(QPainter.Antialiasing)
    left, right, top, bottom = 52, 8, 8, 22
    W, H = pm.width() - left - right, pm.height() - top - bottom
    y0, y1 = float(y.min()), float(y.max())
    if y1 <= y0:
        y1 = y0 + 1
    t0, t1 = float(t[0]), float(t[-1]) if t[-1] > t[0] else float(t[0]) + 1
    p.setPen(QPen(QColor("#999999")))
    p.drawRect(left, top, W, H)
    p.setPen(QPen(QColor("#444444")))
    small = QFont()
    small.setPointSize(8)
    p.setFont(small)
    p.drawText(2, top + 10, f"{y1:.1f}")
    p.drawText(2, top + H, f"{y0:.1f}")
    p.drawText(left, pm.height() - 6, f"{t0:.2f} s")
    p.drawText(left + W - 50, pm.height() - 6, f"{t1:.2f} s")
    p.drawText(left + W // 2 - 60, pm.height() - 6, "heart diameter (µm)")
    pts = QPolygonF([QPointF(left + (ti - t0) / (t1 - t0) * W, top + H - (yi - y0) / (y1 - y0) * H)
                     for ti, yi in zip(t, y)])
    p.setPen(QPen(QColor("#b22222"), 1.2))
    p.drawPolyline(pts)
    p.end()
    return pm


class ReviewPanel(QWidget):
    def __init__(self, get_output_dir: Callable[[], str], run_analysis: Callable[[], None],
                 log: Callable[[str], None]):
        super().__init__()
        self.get_output_dir = get_output_dir
        self.run_analysis = run_analysis
        self.log = log
        self.table_df = None
        self._updating = False
        v = QVBoxLayout(self)

        top = QHBoxLayout()
        refresh = QPushButton("Load traces")
        refresh.setToolTip("Read the traced kymographs and the quality control of the output folder")
        refresh.clicked.connect(self.refresh)
        self.show_combo = QComboBox()
        self.show_combo.addItems(list(SHOW))
        self.show_combo.currentTextChanged.connect(lambda *_: self.fill())
        self.movie_combo = QComboBox()
        self.movie_combo.currentTextChanged.connect(lambda *_: self.fill())
        self.count_lbl = QLabel("")
        top.addWidget(refresh)
        top.addWidget(QLabel("Show:"))
        top.addWidget(self.show_combo)
        top.addWidget(QLabel("Movie:"))
        top.addWidget(self.movie_combo, 1)
        top.addWidget(self.count_lbl)
        v.addLayout(top)

        self.stale = QWidget()
        sl = QHBoxLayout(self.stale)
        sl.setContentsMargins(6, 2, 6, 2)
        self.stale.setStyleSheet("background: #fff3c4; border-radius: 4px")
        sl.addWidget(QLabel("⚠ The selection changed since the last beat analysis - re-run step 3 to update "
                            "the summary tables."), 1)
        rerun = QPushButton("Re-run beat analysis (step 3)")
        rerun.clicked.connect(self.run_analysis)
        sl.addWidget(rerun)
        self.stale.setVisible(False)
        v.addWidget(self.stale)

        split = QSplitter(Qt.Vertical)
        self.table = QTableWidget(0, len(COLS))
        self.table.setHorizontalHeaderLabels(COLS)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(C_NAME, QHeaderView.Stretch)
        self.table.setSortingEnabled(False)
        self.table.itemChanged.connect(self._checkbox_changed)
        self.table.currentCellChanged.connect(lambda *_: self.preview())
        split.addWidget(self.table)

        pv = QWidget()
        pl = QVBoxLayout(pv)
        pl.setContentsMargins(0, 0, 0, 0)
        btns = QHBoxLayout()
        self.add_btn = QPushButton("Add to analysis  (A)")
        self.add_btn.setToolTip("Copy the selected traces into 'excellent traces'")
        self.rm_btn = QPushButton("Remove from analysis  (R)")
        self.rm_btn.setToolTip("Remove the selected traces from 'excellent traces'")
        self.auto_btn = QPushButton("Back to automatic  (U)")
        self.auto_btn.setToolTip("Forget your decision; use the automatic quality control again")
        self.add_btn.clicked.connect(lambda: self.decide("include"))
        self.rm_btn.clicked.connect(lambda: self.decide("exclude"))
        self.auto_btn.clicked.connect(lambda: self.decide("auto"))
        self.stretch = QComboBox()
        self.stretch.addItems(["1×", "2×", "3×", "4×", "6×"])
        self.stretch.setCurrentText("2×")
        self.stretch.setToolTip("Vertical magnification of the kymograph")
        self.stretch.currentTextChanged.connect(lambda *_: self.preview())
        for b in (self.add_btn, self.rm_btn, self.auto_btn):
            btns.addWidget(b)
        btns.addStretch(1)
        btns.addWidget(QLabel("Height:"))
        btns.addWidget(self.stretch)
        pl.addLayout(btns)
        self.info = QLabel("Load the traces, then select one to see the traced kymograph.  "
                           "Keys: A = add, R = remove, U = automatic, ↑/↓ = next trace.")
        self.info.setWordWrap(True)
        pl.addWidget(self.info)
        self.image = QLabel()
        self.image.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        scroll = QScrollArea()
        scroll.setWidget(self.image)
        scroll.setWidgetResizable(True)
        pl.addWidget(scroll, 3)
        self.plot = QLabel()
        self.plot.setMinimumHeight(120)
        pl.addWidget(self.plot, 1)
        split.addWidget(pv)
        split.setSizes([300, 420])
        v.addWidget(split, 1)
        hint = QLabel("<small>* rescued: good trace of a movie without any excellent trace - the automatic QC "
                      "already uses it. Decisions are saved in <i>balled/manual_curation.csv</i> and kept when "
                      "step 2 runs again.</small>")
        hint.setWordWrap(True)
        v.addWidget(hint)

        for key, dec in (("A", "include"), ("R", "exclude"), ("U", "auto")):
            sc = QShortcut(QKeySequence(key), self)
            sc.setContext(Qt.WidgetWithChildrenShortcut)
            sc.activated.connect(lambda d=dec: self.decide(d))

    # ------------------------------------------------------------ data
    def balled(self) -> str:
        return os.path.join(self.get_output_dir(), "balled")

    def refresh(self):
        b = self.balled()
        if not os.path.exists(os.path.join(b, "Quality_control.csv")):
            QMessageBox.information(self, "Review", "No quality control found in\n" + b +
                                    "\n\nRun step 2 (background, tracing and QC) first.")
            return
        self.table_df = curation.trace_table(b)
        movies = ["All movies"] + sorted(self.table_df["movie"].unique())
        cur = self.movie_combo.currentText()
        self.movie_combo.blockSignals(True)
        self.movie_combo.clear()
        self.movie_combo.addItems(movies)
        if cur in movies:
            self.movie_combo.setCurrentText(cur)
        self.movie_combo.blockSignals(False)
        self.fill()

    def visible(self):
        t = self.table_df
        if t is None:
            return t
        sel = SHOW[self.show_combo.currentText()](t)
        m = self.movie_combo.currentText()
        if m and m != "All movies":
            sel &= t["movie"] == m
        return t[sel]

    def fill(self, keep_row: int | None = None):
        t = self.visible()
        if t is None:
            return
        self._updating = True
        self.table.setRowCount(len(t))
        for r, (idx, row) in enumerate(t.iterrows()):
            chk = QTableWidgetItem("")
            chk.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled | Qt.ItemIsSelectable)
            chk.setCheckState(Qt.Checked if row["in_analysis"] else Qt.Unchecked)
            chk.setData(Qt.UserRole, row["csv"])
            self.table.setItem(r, C_IN, chk)
            dec = {"include": "added by you", "exclude": "removed by you"}.get(row["decision"], "–")
            vals = {C_NAME: row["csv"].replace(".tiff.csv", ""), C_MOVIE: row["movie"],
                    C_X: f"{row['Xpos']:.0f}" if row["Xpos"] == row["Xpos"] else "",
                    C_AUTO: AUTO_TEXT.get(row["automatic"], str(row["automatic"])), C_DEC: dec,
                    C_SD: _num(row["sd"]), C_AU: _num(row["autocorr_up"]), C_AD: _num(row["autocorr_down"]),
                    C_DC: _num(row["dist_cor"])}
            for c, text in vals.items():
                it = QTableWidgetItem(text)
                it.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
                self.table.setItem(r, c, it)
            color = "#e3f4e3" if row["in_analysis"] else "#ffffff"
            for c in range(len(COLS)):
                self.table.item(r, c).setBackground(QBrush(QColor(color)))
            if row["decision"]:
                f = QFont()
                f.setBold(True)
                self.table.item(r, C_DEC).setFont(f)
                self.table.item(r, C_DEC).setForeground(QBrush(QColor("#1a4fb3")))
        self._updating = False
        full = self.table_df
        self.count_lbl.setText(f"{len(t)} shown · {int(full['in_analysis'].sum())} of {len(full)} traces in the "
                               f"analysis · {int((full['decision'] != '').sum())} decided by you")
        self.stale.setVisible(curation.analysis_is_stale(self.balled()))
        if len(t):
            self.table.selectRow(min(keep_row if keep_row is not None else 0, len(t) - 1))
            self.preview()     # also when the row number did not change (no currentCellChanged signal)
        else:
            self.image.clear()
            self.plot.clear()

    def selected_csvs(self) -> list[str]:
        rows = sorted({i.row() for i in self.table.selectedIndexes()})
        return [self.table.item(r, C_IN).data(Qt.UserRole) for r in rows if self.table.item(r, C_IN)]

    # ------------------------------------------------------------ actions
    def decide(self, decision: str):
        csvs = self.selected_csvs()
        if not csvs or self.table_df is None:
            return
        row = self.table.currentRow()
        curation.set_decision(self.balled(), csvs, decision)
        word = {"include": "added to", "exclude": "removed from", "auto": "reset to automatic for"}[decision]
        self.log(f"Review: {len(csvs)} trace(s) {word} the analysis")
        self.table_df = curation.trace_table(self.balled())
        # after a single decision go on to the next trace (fast review): if the trace is still
        # listed move one row down, otherwise the next trace has moved up into this row
        still_listed = len(csvs) == 1 and csvs[0] in set(self.visible()["csv"])
        self.fill(keep_row=row + 1 if still_listed else row)

    def _checkbox_changed(self, item):
        if self._updating or item.column() != C_IN:
            return
        csv = item.data(Qt.UserRole)
        want = item.checkState() == Qt.Checked
        self.table.selectRow(item.row())
        curation.set_decision(self.balled(), [csv], "include" if want else "exclude")
        self.log(f"Review: {csv} {'added to' if want else 'removed from'} the analysis")
        self.table_df = curation.trace_table(self.balled())
        self.fill(keep_row=item.row())

    # ------------------------------------------------------------ preview
    def preview(self):
        r = self.table.currentRow()
        it = self.table.item(r, C_IN) if r >= 0 else None
        if it is None or self.table_df is None:
            return
        csv = it.data(Qt.UserRole)
        b = self.balled()
        jpg = os.path.join(b, csv[: -len(".csv")] + "_traced.jpg")
        row = self.table_df[self.table_df["csv"] == csv].iloc[0]
        pm = QPixmap(jpg)
        if not pm.isNull():
            k = int(self.stretch.currentText().rstrip("×"))
            width = max(self.image.parentWidget().width() - 20, 400)
            self.image.setPixmap(pm.scaled(width, pm.height() * k * width // max(pm.width(), 1),
                                           Qt.IgnoreAspectRatio, Qt.SmoothTransformation))
        else:
            self.image.setText("traced image not found: " + os.path.basename(jpg))
        self.plot.setPixmap(diameter_plot(os.path.join(b, csv), max(self.plot.width(), 400),
                                          max(self.plot.height(), 120)))
        state = "IN the analysis" if row["in_analysis"] else "not in the analysis"
        self.info.setText(f"<b>{csv.replace('.tiff.csv', '')}</b> - automatic QC: "
                          f"{AUTO_TEXT.get(row['automatic'], row['automatic'])}, currently <b>{state}</b>. "
                          "Green edges = traced heart walls.  Keys: A add, R remove, U automatic.")


def _num(v) -> str:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return ""
    return "" if f != f else f"{f:.3g}"
