# tdtK Heart Analyzer

A desktop app that analyses movies of fluorescent (tdTomato) _Drosophila melanogaster_ hearts:
it makes kymographs (M-modes), traces the heart walls, and calculates heart period, arrhythmia
index, diastolic/systolic diameters, fractional shortening, contraction velocities and direction
for every fly, summarised per genotype.

It runs on Windows, macOS and Linux and needs nothing else installed: no R, Fiji/ImageJ, Java
or Bio-Formats tools.

[![DOI](https://zenodo.org/badge/192813672.svg)](https://zenodo.org/badge/latestdoi/192813672)

> **Status: 1.0.0 alpha.** Validated on synthetic movies and against reference values; please
> check a few of your own recordings before analysing a whole experiment.

![Analysis tab](tdtk-analyzer/docs/screenshot-analysis.png)

## Features

- **Reads most microscope formats:** Hamamatsu `.cxd`, OME-TIFF, ImageJ/Fiji and Micro-Manager
  TIFF, Zeiss `.czi`/`.lsm`, Nikon `.nd2`, Leica `.lif`, Olympus `.oib`/`.oif`, MetaMorph `.stk`,
  plain TIFF stacks, and video (`.avi`, `.mp4`, `.mov`, `.mkv`).
- **Checks the metadata on import:** missing or doubtful frame intervals and pixel sizes are
  flagged and can be corrected in the app (per movie or for many movies at once).
- **Three steps:** (1) import movies and make kymographs, (2) background subtraction, tracing
  and automatic quality control, (3) beat analysis and summary tables.
- **Manual review:** traces rejected by the quality control can be added back to the analysis
  (and accepted ones removed) with a preview of each trace.
- **Fast:** runs on all CPU cores; interrupted runs continue where they stopped.
- **Command line** for scripted or server use.

## Installation

### Option 1: download the ready-made app (no Python needed)

1. On GitHub, open the **Actions** tab and select the **tdtK Analyzer (Python app)** workflow.
2. Open the most recent run with a green tick and scroll down to **Artifacts**
   (you need to be signed in to GitHub to download them).
3. Download the file for your system:

   | System | Artifact | Contains |
   |---|---|---|
   | Windows 10/11 (64-bit) | `tdtK-Heart-Analyzer-Windows` | folder `tdtK-Heart-Analyzer` |
   | macOS (Apple silicon) | `tdtK-Heart-Analyzer-macOS` | `tdtK-Heart-Analyzer-macOS.tar.gz` |
   | Linux (x86-64) | `tdtK-Heart-Analyzer-Linux` | `tdtK-Heart-Analyzer-Linux.tar.gz` |

4. Unpack and start it:

   - **Windows:** unzip, open the `tdtK-Heart-Analyzer` folder and double-click
     `tdtK-Heart-Analyzer.exe`. Keep the `.exe` together with the `_internal` folder next to it.
     If Windows SmartScreen warns about an unknown publisher, click *More info* → *Run anyway*.
   - **macOS:** unzip, then double-click `tdtK-Heart-Analyzer-macOS.tar.gz` to get
     `tdtK Heart Analyzer.app`, and move it to *Applications*. The app is not signed by Apple, so
     the first time right-click it → *Open* → *Open* (or *System Settings → Privacy & Security →
     Open Anyway*).
   - **Linux:** unzip, then
     ```bash
     tar -xzf tdtK-Heart-Analyzer-Linux.tar.gz
     ./tdtK-Heart-Analyzer/tdtK-Heart-Analyzer
     ```
     A desktop with the usual Qt libraries is needed (on Ubuntu/Debian:
     `sudo apt install libegl1 libgl1 libxkbcommon0 libfontconfig1 libdbus-1-3`).

### Option 2: install with Python

Needs Python 3.9 or newer (tested with 3.9 to 3.13; 3.11 or 3.12 recommended).

```bash
git clone https://github.com/gvogler/FlyHearts-tdtK-Rscripts.git
cd FlyHearts-tdtK-Rscripts/tdtk-analyzer
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install .
tdtk-analyzer-gui                  # start the app
```

With [uv](https://docs.astral.sh/uv/) you can use a newer Python without touching the system one:

```bash
cd FlyHearts-tdtK-Rscripts/tdtk-analyzer
uv venv -p 3.12 .venv && source .venv/bin/activate
uv pip install .
tdtk-analyzer-gui
```

Afterwards, activate the environment (`source .venv/bin/activate`) and run `tdtk-analyzer-gui`
whenever you want to use the app. To update, `git pull` and run `pip install .` again.

Optional: formats not listed above (`.ims`, `.dv`, `.vsi`, `.ics`, `.oir`, …) are read through
Bio-Formats if its command line tools ([bftools](https://www.openmicroscopy.org/bio-formats/downloads/),
needs Java) are on the PATH.

### Option 3: build the app yourself

```bash
cd FlyHearts-tdtK-Rscripts/tdtk-analyzer
pip install -e ".[dev]"
pytest -q                          # optional: run the tests
cd packaging && pyinstaller tdtk-analyzer.spec --noconfirm
```

The app is written to `tdtk-analyzer/packaging/dist/`. PyInstaller builds for the system it runs
on, so build on Windows to get the Windows app.

## Quick start

1. Start **tdtK Heart Analyzer**.
2. On the **Analysis** tab choose
   - the **movie folder** (sub-folders are included),
   - an **output folder**,
   - the **genotype mappings** file (`mappings.xlsx` or `mappings.csv`; an example is in this
     repository).
3. Recommended for formats other than `.cxd`: on the **Import** tab click *Scan movie folder*,
   check the status of each movie and fill in missing frame intervals or pixel sizes.
4. Click **Run analysis**. Progress, log and elapsed time are shown live.
5. Optional: on **Review traces**, add usable M-modes that the quality control rejected, then
   click *Re-run beat analysis (step 3)*.
6. The **Results** tab shows the timing and all summary tables (`final_all_data.csv`,
   `master_table.csv`, `EDD/ESD/FS/HR_table.csv`, …), which are saved in
   `<output folder>/balled/excellent traces`.

Command line equivalent:

```bash
tdtk-analyzer scan --movies /data/movies --output /data/analysis
tdtk-analyzer run  --movies /data/movies --output /data/analysis --mappings /data/mappings.xlsx
```

The full manual (movie formats, metadata checks, review, settings, output files) is in
[tdtk-analyzer/README.md](tdtk-analyzer/README.md).

## License

GPL-3.0, see [LICENSE](LICENSE). The original R scripts are still in this repository; see
[README-R.md](README-R.md).
