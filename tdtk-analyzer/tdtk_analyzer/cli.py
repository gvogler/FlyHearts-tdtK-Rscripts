"""Command line interface: ``tdtk-analyzer run --movies DIR --output DIR --mappings FILE``."""

from __future__ import annotations

import argparse
import sys

from . import __version__
from .pipeline import Callbacks, Pipeline, Settings


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="tdtk-analyzer", description="tdtK fly heart analysis (no R needed).")
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="cmd")
    r = sub.add_parser("run", help="run the analysis")
    r.add_argument("--movies", default="", help="folder with the .cxd movies (searched recursively)")
    r.add_argument("--output", required=True, help="output folder (R: 'target_dir')")
    r.add_argument("--mappings", default="", help="mappings.xlsx or mappings.csv")
    r.add_argument("--steps", default="1,2,3", help="steps to run, e.g. 1,2,3 or 3")
    d = Settings()
    r.add_argument("--workers", type=int, default=d.workers, help="parallel workers for tracing/analysis")
    r.add_argument("--movie-workers", type=int, default=d.movie_workers, help="movies processed in parallel")
    r.add_argument("--min-size-mb", type=float, default=d.min_file_size_mb, help="skip smaller .cxd files")
    r.add_argument("--min-frames", type=int, default=d.min_frames)
    r.add_argument("--max-interval-ms", type=float, default=d.max_frame_interval_ms,
                   help="skip movies with a longer frame interval")
    r.add_argument("--rolling-ball", type=float, default=d.rolling_ball_radius, help="background radius (px)")
    r.add_argument("--channel", default=d.channel, help="'auto', channel index or part of a channel name")
    r.add_argument("--z-plane", default=d.z_plane, help="z plane index or 'max'")
    r.add_argument("--rotate", default=d.rotate, choices=["0", "90", "180", "270", "auto"])
    r.add_argument("--interval-ms", type=float, default=0.0, help="frame interval for files without one")
    r.add_argument("--pixel-um", type=float, default=0.0, help="pixel size for files without one")
    sc = sub.add_parser("scan", help="list the movies and write the editable import table (movie_import.csv)")
    sc.add_argument("--movies", required=True)
    sc.add_argument("--output", required=True, help="folder for movie_import.csv (the analysis output folder)")
    sub.add_parser("gui", help="start the graphical interface")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if args.cmd in (None, "gui"):
        from .gui import main as gui_main
        return gui_main()
    if args.cmd == "scan":
        from .movies import ImportOptions, movie_status
        p = Pipeline(Settings(movie_dir=args.movies, output_dir=args.output))
        infos = p.import_movies()
        for m in infos:
            level, msg = movie_status(m, ImportOptions())
            print(f"{level:5s} {m.rel}{' [' + m.series_name + ']' if m.n_series > 1 else ''}  ({m.format}, "
                  f"{m.size_t}x{m.size_y}x{m.size_x}): {msg}")
        print(f"\n{len(infos)} movies. Edit {args.output}/movie_import.csv to change names, channels, rotation,"
              " frame interval or pixel size; 'run' uses it.")
        return 0
    steps = {s.strip() for s in args.steps.split(",")}
    s = Settings(movie_dir=args.movies, output_dir=args.output, mappings_file=args.mappings,
                 run_kymographs="1" in steps, run_tracing="2" in steps, run_analysis="3" in steps,
                 workers=args.workers, movie_workers=args.movie_workers, min_file_size_mb=args.min_size_mb,
                 min_frames=args.min_frames, max_frame_interval_ms=args.max_interval_ms,
                 rolling_ball_radius=args.rolling_ball, channel=args.channel, z_plane=args.z_plane,
                 rotate=args.rotate, default_interval_ms=args.interval_ms, default_pixel_um=args.pixel_um)
    last = {}

    def progress(stage, done, total):
        if total and (done == total or done - last.get(stage, -1) >= max(1, total // 20)):
            last[stage] = done
            print(f"  {stage}: {done}/{total}", file=sys.stderr)

    try:
        report = Pipeline(s, Callbacks(log=print, progress=progress)).run()
    except ValueError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    except RuntimeError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    print("\nTiming:")
    for t in report.timing:
        print(f"  {t.step:40s} {t.seconds:9.1f} s  ({t.items} items)")
    if report.errors:
        print(f"\n{len(report.errors)} file(s) had errors - see tdtk_analyzer_log.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
