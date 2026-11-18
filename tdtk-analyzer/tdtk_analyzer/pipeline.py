"""Runs the three analysis steps with progress, cancellation, logging and timing.

Used by both the GUI and the command line. Folder layout is the same as the R
script's:

    <movie folder>/**/*.cxd                 -> kymographs, *_directionmarks.csv, metadata
    <output folder>/TIFFs                   copies of the kymographs
    <output folder>/balled                  background-subtracted kymographs + traces
    <output folder>/balled/{excellent,good,bad} traces
    <output folder>/balled/excellent traces summary tables
"""

from __future__ import annotations

import multiprocessing
import os
import shutil
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass, field
from typing import Callable, Iterable

import pandas as pd

from . import __version__, rio
from .aggregate import run_analysis
from .kymograph import KymographSettings, process_movie
from .movies import (ImportOptions, MovieInfo, apply_table, load_table, movie_status, save_table,
                     scan_movies)
from .tracing import quality_control, subtract_background_file, trace_kymograph


@dataclass
class Settings:
    movie_dir: str = ""
    output_dir: str = ""
    mappings_file: str = ""
    run_kymographs: bool = True
    run_tracing: bool = True
    run_analysis: bool = True
    workers: int = max(1, (os.cpu_count() or 2) - 1)
    movie_workers: int = 1              # movies need a lot of memory (whole movie in RAM)
    min_file_size_mb: float = 150.0
    min_frames: int = 200
    max_frame_interval_ms: float = 10.0
    rolling_ball_radius: float = 50.0
    # import (all movie formats)
    channel: str = "auto"             # "auto" (tdTomato-like name, else first), index or name
    z_plane: str = "0"                # index or "max"
    rotate: str = "0"                 # 0/90/180/270 or "auto" (heart must lie along X)
    default_interval_ms: float = 0.0  # used when a file has no frame interval
    default_pixel_um: float = 0.0     # used when a file has no pixel size

    def import_options(self) -> ImportOptions:
        return ImportOptions(channel=self.channel, z_plane=self.z_plane, rotate=self.rotate,
                             default_interval_ms=self.default_interval_ms, default_pixel_um=self.default_pixel_um)

    def kymograph_settings(self) -> KymographSettings:
        return KymographSettings(min_file_size=int(self.min_file_size_mb * 1e6), min_frames=self.min_frames,
                                 max_frame_interval=self.max_frame_interval_ms / 1000.0)


class Cancelled(Exception):
    pass


@dataclass
class Callbacks:
    log: Callable[[str], None] = print
    progress: Callable[[str, int, int], None] = lambda stage, done, total: None
    cancelled: Callable[[], bool] = lambda: False


@dataclass
class TimingRow:
    step: str
    items: int
    seconds: float


@dataclass
class RunReport:
    timing: list = field(default_factory=list)
    counts: dict = field(default_factory=dict)
    errors: list = field(default_factory=list)


def _parallel(func, tasks: list, workers: int, stage: str, cb: Callbacks) -> list:
    """Run func(task) for all tasks; results in task order; reports progress."""
    total = len(tasks)
    cb.progress(stage, 0, total)
    if total == 0:
        return []
    results = [None] * total
    if workers <= 1 or total == 1:
        for i, t in enumerate(tasks):
            if cb.cancelled():
                raise Cancelled()
            results[i] = func(t)
            cb.progress(stage, i + 1, total)
        return results
    # 'spawn' everywhere: fork() of a multi-threaded process (the GUI runs the pipeline in a
    # thread) can deadlock on Linux; Windows and macOS use spawn already
    with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("spawn")) as ex:
        futs = {ex.submit(func, t): i for i, t in enumerate(tasks)}
        done = 0
        try:
            for f in as_completed(futs):
                results[futs[f]] = f.result()
                done += 1
                cb.progress(stage, done, total)
                if cb.cancelled():
                    raise Cancelled()
        except Cancelled:
            for f in futs:
                f.cancel()
            raise
    return results


def _star(args):
    func, a = args
    return func(*a)


class Pipeline:
    def __init__(self, settings: Settings, callbacks: Callbacks | None = None):
        self.s = settings
        self.cb = callbacks or Callbacks()
        self.report = RunReport()
        self._logfile = None

    # ------------------------------------------------------------------ helpers
    def log(self, msg: str) -> None:
        line = time.strftime("%H:%M:%S ") + msg
        self.cb.log(line)
        if self._logfile:
            self._logfile.write(line + "\n")
            self._logfile.flush()

    def _timed(self, name: str, fn: Callable[[], int]) -> None:
        self.log(f"=== {name} ===")
        t0 = time.perf_counter()
        n = fn()
        dt = time.perf_counter() - t0
        self.report.timing.append(TimingRow(name, n, dt))
        per = f", {dt / n:.2f} s per item" if n else ""
        self.log(f"--- {name}: {dt:.1f} s ({n} items{per})")

    def _results(self, results: Iterable, what: str) -> None:
        for r in results:
            if r is None:
                continue
            if r.status == "error":
                self.report.errors.append((what, r.name, r.message))
                self.log(f"  ERROR {r.name}: {r.message}")
            elif r.status == "skipped" and r.message not in ("already processed", "already traced"):
                self.log(f"  skipped {r.name}: {r.message}")

    # ------------------------------------------------------------------ steps
    def import_movies(self) -> list[MovieInfo]:
        """Find and describe all movies; apply and refresh the import table (movie_import.csv)."""
        def prog(done, total, name):
            self.cb.progress("Import", done, total)
            if self.cb.cancelled():
                raise Cancelled()

        infos = scan_movies(self.s.movie_dir, self.s.output_dir, prog)
        table = load_table(self.s.output_dir) if self.s.output_dir else None
        if table is not None:
            apply_table(infos, table)
            self.log(f"Import table applied: {os.path.join(self.s.output_dir, 'movie_import.csv')}")
        if self.s.output_dir:
            save_table(infos, self.s.output_dir)
        return infos

    def step_kymographs(self) -> int:
        infos = self.import_movies()
        opts, ks = self.s.import_options(), self.s.kymograph_settings()
        formats = pd.Series([m.format for m in infos]).value_counts().to_dict() if infos else {}
        self.log(f"{len(infos)} movies found in {self.s.movie_dir}"
                 + (" (" + ", ".join(f"{k}: {v}" for k, v in formats.items()) + ")" if formats else ""))
        todo = []
        for m in infos:
            level, msg = movie_status(m, opts, ks.min_frames, ks.max_frame_interval)
            if level in ("error", "skip") and msg != "not selected":
                self.log(f"  {'ERROR' if level == 'error' else 'skipped'} {m.rel} [{m.name}]: {msg}")
                if level == "error":
                    self.report.errors.append(("import", m.rel, msg))
            elif level != "skip":
                todo.append(m)
        tasks = [(process_movie, (m, opts, ks)) for m in todo]
        res = _parallel(_star, tasks, self.s.movie_workers, "Kymographs", self.cb)
        self._results(res, "kymographs")
        self.report.counts["movies"] = sum(1 for r in res if r.status == "done")
        return len(infos)

    def step_tracing(self) -> int:
        out = self.s.output_dir
        tiffs_dir, balled = os.path.join(out, "TIFFs"), os.path.join(out, "balled")
        os.makedirs(tiffs_dir, exist_ok=True)
        os.makedirs(balled, exist_ok=True)
        csvs = rio.list_files(self.s.movie_dir, r"csv$", recursive=True, full_names=True)
        if not csvs:
            raise RuntimeError("No CSV files found. Please select a different folder")
        for f in csvs:
            dst = os.path.join(balled, os.path.basename(f))
            if not os.path.exists(dst):
                shutil.copy2(f, dst)
        kymos = rio.list_files(self.s.movie_dir, r"Xpos", recursive=True, full_names=True)
        if not kymos:
            raise RuntimeError("No TIFF files found. Please select a different folder")
        for f in kymos:
            dst = os.path.join(tiffs_dir, os.path.basename(f))
            if not os.path.exists(dst):
                shutil.copy2(f, dst)

        names = rio.list_files(tiffs_dir, r"\.tiff$")
        todo = [n for n in names if not os.path.exists(os.path.join(balled, n))]
        self.log(f"Background subtraction (rolling ball {self.s.rolling_ball_radius:g} px) on {len(todo)} kymographs")
        tasks = [(subtract_background_file, (os.path.join(tiffs_dir, n), os.path.join(balled, n),
                                             self.s.rolling_ball_radius)) for n in todo]
        self._results(_parallel(_star, tasks, self.s.workers, "Background", self.cb), "background")

        names = rio.list_files(balled, r"\.tiff$")
        self.log(f"Tracing {len(names)} kymographs")
        tasks = [(trace_kymograph, (n, balled)) for n in names]
        self._results(_parallel(_star, tasks, self.s.workers, "Tracing", self.cb), "tracing")

        self.log("Quality control")
        qc = quality_control(balled)
        if not qc.empty:
            counts = qc["quality"].value_counts().to_dict()
            self.report.counts["quality"] = counts
            self.log("  " + ", ".join(f"{k}: {v}" for k, v in counts.items()))
        return len(names)

    def step_analysis(self) -> int:
        excellent = os.path.join(self.s.output_dir, "balled", "excellent traces")
        if not os.path.isdir(excellent):
            raise RuntimeError(f"Folder not found: {excellent} (run the tracing step first)")

        def map_fn(func, tasks):
            tasks = list(tasks)
            return _parallel(func, tasks, self.s.workers, "Beat analysis", self.cb)

        summary = run_analysis(excellent, self.s.mappings_file, map_fn=map_fn, log=self.log)
        self.report.counts["analysis"] = {k: v for k, v in summary.items() if k != "skipped"}
        for name, msg in summary["skipped"]:
            self.report.errors.append(("analysis", name, msg))
        self.log(f"  {summary['analyzed']} of {summary['kymographs']} kymographs analyzed, "
                 f"{summary['genotypes']} genotype groups")
        return summary["kymographs"]

    # ------------------------------------------------------------------ run
    def validate(self) -> list[str]:
        problems = []
        if (self.s.run_kymographs or self.s.run_tracing) and not os.path.isdir(self.s.movie_dir):
            problems.append("Select the folder with the movies.")
        if not self.s.output_dir:
            problems.append("Select an output folder.")
        if self.s.run_analysis:
            name = os.path.basename(self.s.mappings_file).lower()
            if not os.path.isfile(self.s.mappings_file):
                problems.append("Select the genotype mappings file (mappings.xlsx or mappings.csv).")
            elif not name.endswith((".xlsx", ".xls", ".csv")):
                problems.append("The mappings file must be an .xlsx or .csv file.")
        if not (self.s.run_kymographs or self.s.run_tracing or self.s.run_analysis):
            problems.append("Select at least one step.")
        return problems

    def run(self) -> RunReport:
        problems = self.validate()
        if problems:
            raise ValueError("\n".join(problems))
        os.makedirs(self.s.output_dir, exist_ok=True)
        self._logfile = open(os.path.join(self.s.output_dir, "tdtk_analyzer_log.txt"), "a", encoding="utf-8")
        try:
            self.log(f"tdtK Analyzer {__version__}")
            for k, v in asdict(self.s).items():
                self.log(f"  {k:22s} {v}")
            t0 = time.perf_counter()
            if self.s.run_kymographs:
                self._timed("1 Import movies and make kymographs", self.step_kymographs)
            if self.s.run_tracing:
                self._timed("2 Background, tracing and QC", self.step_tracing)
            if self.s.run_analysis:
                self._timed("3 Beat analysis and summary tables", self.step_analysis)
            total = time.perf_counter() - t0
            self.report.timing.append(TimingRow("Total", 0, total))
            pd.DataFrame([asdict(t) for t in self.report.timing]).to_csv(
                os.path.join(self.s.output_dir, "timing.csv"), index=False)
            self.log(f"Finished in {total:.1f} s. Timing saved to timing.csv")
            return self.report
        except Cancelled:
            self.log("Cancelled.")
            raise
        finally:
            self._logfile.close()
            self._logfile = None
