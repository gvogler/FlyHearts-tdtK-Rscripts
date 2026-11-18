"""Manual review: adding 'good'/'bad' traces to the analysis and removing excellent ones."""

import os

import pandas as pd

from synthetic import heart_movie, write_cxd
from tdtk_analyzer import curation
from tdtk_analyzer.pipeline import Pipeline, Settings
from tdtk_analyzer.tracing import quality_control

MAPPINGS = ("CODE,cross,type,Gal4-line,UAS-line,Gene,human ortholog,DIOPT,Stock collection\n"
            "MAYO0001,Hand x ctrl,control,Hand,ctrl,CG1,GENE1,8 of 11,BL\n"
            "MAYO0002,Hand x rnai,experiment,Hand,rnai,CG2,GENE2,8 of 11,VDRC\n")


def _run(tmp_path):
    movies = tmp_path / "movies"
    movies.mkdir()
    for name, seed, rev, per in [("MAYO0001_1_1wf_a.cxd", 0, False, 0.25), ("MAYO0001_2_1wf_a.cxd", 1, False, 0.27),
                                 ("MAYO0002_1_1wf_a.cxd", 2, True, 0.30), ("MAYO0002_2_1wf_a.cxd", 3, True, 0.22)]:
        write_cxd(str(movies / name), heart_movie(seed=seed, reverse=rev, period=per))
    (tmp_path / "mappings.csv").write_text(MAPPINGS)
    s = Settings(movie_dir=str(movies), output_dir=str(tmp_path / "out"),
                 mappings_file=str(tmp_path / "mappings.csv"), workers=1, min_file_size_mb=1)
    Pipeline(s).run()
    return s, str(tmp_path / "out" / "balled")


def test_include_exclude_and_rerun(tmp_path):
    s, balled = _run(tmp_path)
    exc = os.path.join(balled, "excellent traces")
    t = curation.trace_table(balled)
    assert set(t["automatic"]) <= {"excellent", "rescued", "good", "bad"}
    outside = t[~t["in_analysis"]].sort_values("automatic", key=lambda a: a != "good")["csv"].tolist()
    inside = t[t["in_analysis"]]["csv"].tolist()
    assert outside and inside

    # add a rejected trace, remove an accepted one
    add, drop = outside[0], inside[0]
    curation.set_decision(balled, [add], "include")
    curation.set_decision(balled, [drop], "exclude")
    assert os.path.exists(os.path.join(exc, add)) and os.path.exists(os.path.join(exc, add[:-4] + "_traced.jpg"))
    assert not os.path.exists(os.path.join(exc, drop))
    assert curation.analysis_is_stale(balled)
    dec = pd.read_csv(os.path.join(balled, curation.CURATION))
    assert dict(zip(dec["csv"], dec["decision"])) == {add: "include", drop: "exclude"}

    # decisions survive a new quality control run (step 2)
    quality_control(balled)
    assert os.path.exists(os.path.join(exc, add)) and not os.path.exists(os.path.join(exc, drop))

    # step 3 uses the reviewed selection
    s.run_kymographs = s.run_tracing = False
    Pipeline(s).run()
    analyzed = set(pd.read_csv(os.path.join(exc, "all_data_table.csv"))["file"])
    assert add[:-4] + "_traced.jpg" in analyzed           # the added trace is analyzed ...
    assert drop[:-4] + "_traced.jpg" not in analyzed      # ... the removed one is not
    assert not curation.analysis_is_stale(balled)

    # back to automatic restores the QC's choice
    curation.set_decision(balled, [add, drop], "auto")
    assert not os.path.exists(os.path.join(exc, add)) and os.path.exists(os.path.join(exc, drop))
    assert curation.load_decisions(balled) == {}
