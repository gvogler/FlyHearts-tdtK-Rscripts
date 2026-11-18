"""Checks a built app (PyInstaller output) end to end, without Python packages from the build.

    python packaging/smoke_test.py packaging/dist/tdtK-Heart-Analyzer

1. writes synthetic heart movies in five formats (CXD, OME-TIFF, ImageJ TIFF, CZI, plain TIFF)
   with known beat periods, plus an AVI,
2. runs 'tdtk-analyzer-cli scan' and 'run' (all three steps) from the built folder,
3. checks the summary tables (every movie analysed, heart periods correct),
4. starts the GUI executable in self-test mode (opens the main window and quits).
Exits with a non-zero code on the first problem.
"""

import os
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tests"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import tifffile  # noqa: E402
from synthetic import heart_movie, write_cxd, write_imagej_avi  # noqa: E402

MAPPINGS = ("CODE,cross,type,Gal4-line,UAS-line,Gene,human ortholog,DIOPT,Stock collection\n"
            "MAYO0001,Hand x ctrl,control,Hand,ctrl,CG1,GENE1,8 of 11,BL\n"
            "MAYO0002,Hand x rnai,experiment,Hand,rnai,CG2,GENE2,8 of 11,VDRC\n")
PERIOD = {"MAYO0001_1_1wf": 0.25, "MAYO0001_2_1wf": 0.27, "MAYO0002_1_1wf": 0.30,
          "MAYO0002_2_1wf": 0.22, "MAYO0002_3_1wf": 0.26}


def exe(folder: str, name: str) -> str:
    path = os.path.join(folder, name + (".exe" if sys.platform == "win32" else ""))
    if not os.path.isfile(path):
        sys.exit(f"FAIL: {path} not found")
    return path


def write_movies(mdir: str) -> None:
    from pylibCZIrw import czi as pyczi

    mv = {k: heart_movie(seed=i, period=p, reverse=k.startswith("MAYO0002")) for i, (k, p) in enumerate(PERIOD.items())}
    write_cxd(os.path.join(mdir, "MAYO0001_1_1wf_a.cxd"), mv["MAYO0001_1_1wf"])
    both = np.stack([mv["MAYO0001_2_1wf"] // 3, mv["MAYO0001_2_1wf"]], axis=1).astype(np.uint16)
    tifffile.imwrite(os.path.join(mdir, "MAYO0001_2_1wf_a.ome.tif"), both, ome=True, metadata={
        "axes": "TCYX", "PhysicalSizeX": 0.65, "TimeIncrement": 0.005, "Channel": {"Name": ["GFP", "tdTomato"]}})
    tifffile.imwrite(os.path.join(mdir, "MAYO0002_1_1wf_a.tif"), mv["MAYO0002_1_1wf"], imagej=True,
                     resolution=(1 / 0.65, 1 / 0.65), metadata={"axes": "TYX", "finterval": 0.005, "unit": "um"})
    with pyczi.create_czi(os.path.join(mdir, "MAYO0002_2_1wf_a.czi"), exist_ok=True) as w:
        for t, fr in enumerate(mv["MAYO0002_2_1wf"]):
            w.write(data=fr.astype(np.uint16)[..., None], plane={"T": t, "C": 0, "Z": 0})
        w.write_metadata(document_name="x", channel_names={0: "tdTomato"}, scale_x=0.65e-6, scale_y=0.65e-6)
    tifffile.imwrite(os.path.join(mdir, "MAYO0002_3_1wf_a.tif"), mv["MAYO0002_3_1wf"])   # no metadata
    avi_dir = os.path.join(mdir, "..", "avi")
    os.makedirs(avi_dir, exist_ok=True)
    write_imagej_avi(os.path.join(avi_dir, "MAYO0003_1_1wf_a.avi"), mv["MAYO0001_1_1wf"])


def run(cmd: list, timeout: int = 900, env=None) -> str:
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=env)
    print(f"$ {' '.join(os.path.basename(c) if i == 0 else c for i, c in enumerate(cmd))}  "
          f"[exit {p.returncode}, {time.time() - t0:.1f} s]")
    out = (p.stdout or "") + (p.stderr or "")
    if p.returncode != 0:
        print(out[-4000:])
        sys.exit(f"FAIL: exit code {p.returncode}")
    return out


def main() -> int:
    folder = os.path.abspath(sys.argv[1])
    cli, gui = exe(folder, "tdtk-analyzer-cli"), exe(folder, "tdtK-Heart-Analyzer")
    print(run([cli, "--version"]).strip())
    with tempfile.TemporaryDirectory() as tmp:
        mdir, out = os.path.join(tmp, "movies"), os.path.join(tmp, "analysis")
        os.makedirs(mdir)
        write_movies(mdir)
        maps = os.path.join(tmp, "mappings.csv")
        with open(maps, "w") as fh:
            fh.write(MAPPINGS)

        # every movie reader is part of the build (AVI via the video library)
        scan = run([cli, "scan", "--movies", os.path.join(tmp, "avi"), "--output", os.path.join(tmp, "avi_out")])
        print(scan.strip())
        line = next((ln for ln in scan.splitlines() if "MAYO0003_1_1wf_a.avi" in ln), "")
        if "500x56x240" not in line:          # frames x height x width read from the file
            sys.exit("FAIL: the AVI file was not read")

        run([cli, "run", "--movies", mdir, "--output", out, "--mappings", maps, "--min-size-mb", "1",
             "--workers", "2", "--interval-ms", "5", "--pixel-um", "0.65"])
        ex = os.path.join(out, "balled", "excellent traces")
        for f in ("master_table.csv", "final_all_data.csv", "EDD_table.csv", "HR_table.csv"):
            if not os.path.isfile(os.path.join(ex, f)):
                sys.exit(f"FAIL: {f} was not written")
        for f in os.listdir(mdir):
            if f.endswith(("cxd", "tif", "czi")) and not os.path.exists(os.path.join(mdir, f + "_directionmarks.csv")):
                sys.exit(f"FAIL: {f} was not analysed")
        fad = pd.read_csv(os.path.join(ex, "final_all_data.csv"))
        if fad.empty:
            sys.exit("FAIL: final_all_data.csv is empty")
        for _, r in fad.iterrows():
            err = abs(r["HP_median"] - PERIOD[r["flyID"]])
            print(f"  {r['flyID']}: heart period {r['HP_median']:.3f} s (expected {PERIOD[r['flyID']]:.3f})")
            if err >= 0.01:
                sys.exit(f"FAIL: wrong heart period for {r['flyID']}")
        print(f"analysis OK: {len(fad)} traces of {fad['flyID'].nunique()} flies, {fad['cxd_file'].nunique()} movie formats in the summary")

    env = dict(os.environ, TDTK_ANALYZER_SELFTEST="1")
    if sys.platform.startswith("linux"):
        env.setdefault("QT_QPA_PLATFORM", "offscreen")
    run([gui], timeout=180, env=env)
    print("GUI OK: main window opened and closed")
    print("SMOKE TEST PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
