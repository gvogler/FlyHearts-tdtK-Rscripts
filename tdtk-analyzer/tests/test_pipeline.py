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
