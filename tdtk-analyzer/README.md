# tdtK Heart Analyzer (standalone)

A standalone desktop app that runs the same analysis as
`tdtK_Full_Analysis_script_v0.7.5.R`, **without R, RStudio, Fiji/ImageJ, Java or bftools**.
It is a Python port with a graphical interface, and it can be packaged as a double-clickable
app for Windows, macOS and Linux.

> **Status: 1.0.0 alpha.** The port has been validated component by component (see
> *Validation*) and end-to-end on synthetic movies. It has **not yet been compared against
> R results on real `.cxd` recordings**. Please do that before using it for real data (see
> *How to compare with the R script*).

![Analysis tab](docs/screenshot-analysis.png)

## What it does

The same three stages as the R script, with the same folder layout and output files:

| Step | R script part | What happens |
|---|---|---|
| 1 | Script No.1 | Reads every `.cxd` movie (sub-folders included), finds the heart stripes, writes kymograph TIFFs, `…_SD and peaklines.tiff`, `…_new_meta_data.csv` and `…_directionmarks.csv` next to each movie |
| 2 | Fiji macro, Scripts No.2 and No.3 | Copies the kymographs to `<output>/TIFFs`, subtracts the background (ImageJ rolling ball, radius 50), traces the heart edges, writes `.tiff.csv` and `_traced.jpg`, and sorts traces into *excellent*, *good* and *bad traces* (quality control) |
| 3 | Scripts No.4 and No.5 | Finds beats, fits splines, calculates intervals, arrhythmia index, diameters, fractional shortening, velocities and direction, and writes all summary tables into `<output>/balled/excellent traces` |

Files that already exist are not recomputed, as in the R script, so an interrupted run can
simply be started again.

## Using the app

1. Start *tdtK Heart Analyzer* (or `python -m tdtk_analyzer` from this folder).
2. **Movie folder**: the folder with the `.cxd` files (R: `movie_dir`).
   **Output folder**: R's `target_dir`.
   **Genotype mappings**: `mappings.xlsx` or `mappings.csv`, in the same format as before.
3. Tick the steps to run. These replace R's two "reprocess CXD/TIFF files?" questions.
4. Click **Run analysis**. The progress bars, the log and the elapsed time update live.
   **Cancel** stops after the files that are currently being processed.
5. The **Results** tab shows how long each step took, all the summary tables and a preview of
   every traced kymograph. Double-click a file to open it.

The **Settings** tab has the number of parallel workers and the skip rules used by the R script
(minimum movie size 150 MB, at least 200 frames, frame interval of at most 10 ms), plus the
rolling-ball radius. The defaults are the R script's values.

![Results tab](docs/screenshot-results.png)

### Command line

```bash
python -m tdtk_analyzer run --movies /data/movies --output /data/analysis \
       --mappings /data/movies/mappings.xlsx --workers 8
# only redo the summary tables:
python -m tdtk_analyzer run --output /data/analysis --mappings mappings.xlsx --steps 3
```

Every run appends to `tdtk_analyzer_log.txt` and writes `timing.csv` into the output folder.

## Installing

**Pre-built app:** download the artifact for your system from the *tdtK Analyzer (Python app)*
GitHub Actions workflow, unzip it and start `tdtK-Heart-Analyzer`. Python is not needed.

**From source** (Python 3.10 or newer):

```bash
cd tdtk-analyzer
pip install .            # or: pip install -e ".[dev]" for tests and packaging
tdtk-analyzer-gui        # GUI
tdtk-analyzer run ...    # command line
```

**Building the app yourself:** run `pip install -e ".[dev]"`, then
`cd packaging && pyinstaller tdtk-analyzer.spec --noconfirm`. The app is written to
`packaging/dist/`. PyInstaller builds for the system it runs on, so build on Windows to get a
Windows app. The GitHub workflow builds all three systems.

## Speed

- The CXD reader loads a movie as 16-bit integers (R used 8-byte doubles), so a movie needs
  about a quarter of the memory it needed in R.
- A 157 MB synthetic movie (512×128 px, 1200 frames) goes through step 1 in about 2 s.
- Tracing, quality control and beat analysis run in parallel across all CPU cores.
- Movies are processed one at a time by default, because each is held in memory. Increase
  *Movies processed at the same time* only if you have enough RAM.

## Output files

The same names and columns as the R script, including the column order of `final_all_data.csv`,
`final_all_data_per_fly*.csv`, `all_data_table.csv`, `master_table.csv`, `EDD/ESD/FS/HR_table.csv`,
`min_velocity_data.csv`, `max_velocity_table.csv`, `meta_data_all.csv` and `Quality_control.csv`.

The R script also saved intermediate `.Rdata` files. Instead, the app writes:

| R | App |
|---|---|
| `intervals_final.Rdata` | `intervals_final.csv` |
| `final_transients.Rdata` | `final_transients.csv` |
| `transients_final.Rdata` | `transients_final.csv.gz` |
| `metrics_all_splined.Rdata`, `spline_all.Rdata`, `z.Rdata`, other `.Rdata` | not written (their content is in the tables above) |

## How the R dependencies were replaced

| R / external tool | Replacement | Validation |
|---|---|---|
| bftools `showinf`, RBioFormats `read.image` (CXD) | `cxd.py`: reads the OLE2 file directly, following Bio-Formats' `PCIReader` (frame order, padded rows, timestamps, calibration) | Round-trip test with generated CXD files |
| Fiji macro `Subtract Background… rolling=50` | `background.py`: port of ImageJ's `BackgroundSubtracter` (3×3 smoothing, shrink ×4, ball, bilinear enlarge) | Equal to a line-by-line port of ImageJ's Java loops (tests) |
| `stats::smooth.spline` (GCV) | `rstats.SmoothSpline`: R's knot rule, R's penalty matrix (including sgram's `.3330` constant) and R's Brent search for `spar` | Matches real R: `spar` to 1e-8, predictions to about 1e-12 (tests, reference data in `tests/data`) |
| `baseline::rollingBall(wm=20, ws=20)` | line-by-line port | Matches real R to 1e-15 |
| `EBImage::gblur`, `writeImage`, zoo, pracma, `lm`, `quantile`, `mad`, … | NumPy/SciPy equivalents with R's definitions | Unit tests against values produced by R |

## Differences from the R script

- **Problem files are skipped, not fatal.** A file that would have stopped the whole R run with
  an error (for example no calibration in the CXD file, no complete beats, edges that cannot be
  traced, or a missing direction file) is skipped and reported in the log instead.
- The loops that adapt the peak-window size are capped at 1000 rounds. In R the second loop
  could run forever when the window size reached 0.
- Genotype groups are sorted case-insensitively, which approximates R's locale sorting. This
  only affects the row order.
- **R quirks that are kept on purpose, so the numbers stay comparable:**
  - The direction speed uses `Xpos2 - Xpos1 * resolution`, R's operator precedence.
  - The magnification-2 calibration "fix" has no effect.
  - `final_all_data_per_fly_anterior/posterior.csv` contain `anterograde_speed.x/.y`, where the
    main table has `stocks`.

## How to compare with the R script

1. Copy a few `.cxd` movies into two separate folders.
2. Run R v0.7.5 on one copy and this app on the other, with the same `mappings` file.
3. Compare `Quality_control.csv` (were the same traces classified the same way?), then
   `final_all_data.csv` and `master_table.csv`.

Small differences in the last digits are expected: TIFF/JPEG encoding, and float32 against
float64 in the background subtraction. Please report anything larger.

## Code overview

```
tdtk_analyzer/
  cxd.py         CXD reader (no Java)
  kymograph.py   step 1
  background.py  ImageJ rolling-ball background subtraction
  tracing.py     step 2 (tracing + quality control)
  transients.py  step 3, per kymograph (beats, splines, metrics)
  aggregate.py   step 3, summary tables
  rstats.py      R-compatible numerics (smooth.spline, rollingBall, quantile, lm, ...)
  rio.py         R-style CSV/TIFF/JPEG writing, R-style file naming helpers
  pipeline.py    runs the steps: parallel workers, progress, cancel, log, timing
  gui.py / cli.py
tests/           pytest suite (synthetic CXD writer, R reference values)
packaging/       PyInstaller spec
```
