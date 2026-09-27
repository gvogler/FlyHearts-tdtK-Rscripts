"""End-to-end run on synthetic movies with known beat period and direction."""

import os

import pandas as pd

from synthetic import heart_movie, write_cxd
from tdtk_analyzer.pipeline import Pipeline, Settings

MAPPINGS = ("CODE,cross,type,Gal4-line,UAS-line,Gene,human ortholog,DIOPT,Stock collection\n"
            "MAYO0001,Hand x ctrl,control,Hand,ctrl,CG1,GENE1,8 of 11,BL\n"
            "MAYO0002,Hand x rnai,experiment,Hand,rnai,CG2,GENE2,8 of 11,VDRC\n")


def test_full_pipeline(tmp_path):
    movies = tmp_path / "movies" / "day1"
    movies.mkdir(parents=True)
    specs = [("MAYO0001_1_1wf_a.cxd", 0, False, 0.25), ("MAYO0001_2_1wf_a.cxd", 1, False, 0.27),
             ("MAYO0002_1_1wf_a.cxd", 2, True, 0.30), ("MAYO0002_2_1wf_a.cxd", 3, True, 0.22)]
    for name, seed, rev, per in specs:
        write_cxd(str(movies / name), heart_movie(seed=seed, reverse=rev, period=per))
    (tmp_path / "mappings.csv").write_text(MAPPINGS)
    s = Settings(movie_dir=str(tmp_path / "movies"), output_dir=str(tmp_path / "out"),
                 mappings_file=str(tmp_path / "mappings.csv"), workers=2, min_file_size_mb=1)
    report = Pipeline(s).run()
    ex = tmp_path / "out" / "balled" / "excellent traces"
    for f in ("master_table.csv", "final_all_data.csv", "final_all_data_per_fly.csv", "all_data_table.csv",
              "EDD_table.csv", "HR_table.csv", "meta_data_all.csv", "final_transients.csv"):
        assert (ex / f).exists(), f
    fad = pd.read_csv(ex / "final_all_data.csv")
    assert len(fad.columns) == 57 and list(fad.columns[:3]) == ["CODE", "CODE_long", "flyID"]
    true_period = {"MAYO0001_1_1wf": 0.25, "MAYO0001_2_1wf": 0.27, "MAYO0002_1_1wf": 0.30, "MAYO0002_2_1wf": 0.22}
    for _, row in fad.iterrows():
        assert abs(row["HP_median"] - true_period[row["flyID"]]) < 0.01
    fwd = fad[fad["CODE"] == "MAYO0001"]["anterograde_percent"]
    assert (fwd == 0).all()                              # waves travel towards larger X
    assert os.path.exists(tmp_path / "out" / "timing.csv")
    assert [t.step for t in report.timing][-1] == "Total"


def test_mixed_formats(tmp_path):
    """The same analysis on OME-TIFF, ImageJ TIFF, CZI, CXD and a plain TIFF configured via the import table."""
    import numpy as np
    import tifffile
    from pylibCZIrw import czi as pyczi

    from tdtk_analyzer.movies import MANIFEST

    mdir = tmp_path / "movies"
    mdir.mkdir()
    period = {"MAYO0001_1_1wf": 0.25, "MAYO0001_2_1wf": 0.27, "MAYO0002_1_1wf": 0.30,
              "MAYO0002_2_1wf": 0.22, "MAYO0002_3_1wf": 0.26}
    mv = {k: heart_movie(seed=i, period=p, reverse=k.startswith("MAYO0002")) for i, (k, p) in enumerate(period.items())}
    write_cxd(str(mdir / "MAYO0001_1_1wf_a.cxd"), mv["MAYO0001_1_1wf"])
    both = np.stack([mv["MAYO0001_2_1wf"] // 3, mv["MAYO0001_2_1wf"]], axis=1).astype(np.uint16)
    tifffile.imwrite(mdir / "MAYO0001_2_1wf_a.ome.tif", both, ome=True, metadata={
        "axes": "TCYX", "PhysicalSizeX": 0.65, "TimeIncrement": 0.005, "Channel": {"Name": ["GFP", "tdTomato"]}})
    tifffile.imwrite(mdir / "MAYO0002_1_1wf_a.tif", mv["MAYO0002_1_1wf"], imagej=True,
                     resolution=(1 / 0.65, 1 / 0.65), metadata={"axes": "TYX", "finterval": 0.005, "unit": "um"})
    with pyczi.create_czi(str(mdir / "MAYO0002_2_1wf_a.czi"), exist_ok=True) as w:
        for t, fr in enumerate(mv["MAYO0002_2_1wf"]):
            w.write(data=fr.astype(np.uint16)[..., None], plane={"T": t, "C": 0, "Z": 0})
        w.write_metadata(document_name="x", channel_names={0: "tdTomato"}, scale_x=0.65e-6, scale_y=0.65e-6)
    tifffile.imwrite(mdir / "MAYO0002_3_1wf_a.tif", mv["MAYO0002_3_1wf"])          # no metadata at all
    (tmp_path / "mappings.csv").write_text(MAPPINGS)
    out = tmp_path / "out"
    s = Settings(movie_dir=str(mdir), output_dir=str(out), mappings_file=str(tmp_path / "mappings.csv"),
                 workers=2, min_file_size_mb=1)

    # 1) import only: the plain TIFF lacks frame interval / pixel size
    p = Pipeline(s)
    table = __import__("tdtk_analyzer.movies", fromlist=["to_table"]).to_table(p.import_movies())
    assert set(table["format"]) >= {"Hamamatsu CXD", "OME-TIFF", "ImageJ TIFF", "Zeiss CZI", "TIFF stack"}
    row = table["file"] == "MAYO0002_3_1wf_a.tif"
    assert table.loc[row, "frame_interval_ms"].isna().all()
    # the CZI written here has no time increment either -> fill both in the import table
    table.loc[row | (table["format"] == "Zeiss CZI"), "frame_interval_ms"] = 5
    table.loc[row, "pixel_size_um"] = 0.65
    table.to_csv(out / MANIFEST, index=False)

    # 2) full run
    Pipeline(s).run()
    for f in ("MAYO0001_1_1wf_a.cxd", "MAYO0001_2_1wf_a.ome.tif", "MAYO0002_1_1wf_a.tif",
              "MAYO0002_2_1wf_a.czi", "MAYO0002_3_1wf_a.tif"):
        assert (mdir / (f + "_directionmarks.csv")).exists(), f        # every movie was analyzed
    fad = pd.read_csv(out / "balled" / "excellent traces" / "final_all_data.csv")
    assert fad["cxd_file"].nunique() >= 3                                  # several formats passed QC
    for _, r in fad.iterrows():
        assert abs(r["HP_median"] - period[r["flyID"]]) < 0.01, r["flyID"]
    assert set(fad["cxd_file"]) <= {"MAYO0001_1_1wf_a.cxd", "MAYO0001_2_1wf_a.ome.tif", "MAYO0002_1_1wf_a.tif",
                                    "MAYO0002_2_1wf_a.czi", "MAYO0002_3_1wf_a.tif"}
