from tdtk_analyzer import rio
from tdtk_analyzer.aggregate import age_of, comment_of, fly_part, genotype_of, id_of, no_peak_number, sex_of, xpos_of

F = "MAYO0001_12_1wf_a_b_note.cxd_peak_2_at Xpos_130.tiff_traced.jpg"


def test_filename_parsing_like_r():
    assert genotype_of(F) == "MAYO0001"
    assert id_of(F) == "12"
    assert age_of(F) == "1"
    assert sex_of(F) == "f"
    assert fly_part(F) == "12_1wf"
    assert comment_of(F) == "note"
    assert xpos_of(F) == 130
    assert no_peak_number(F) == "MAYO0001_12_1wf_a_b_note.cxd_peak_X_at Xpos_130.tiff_traced.jpg"
    assert rio.cxd_prefix("d/x.cxd") == "d/x."
    assert rio.cxd_stem(F) == "MAYO0001_12_1wf_a_b_note.cxd"


def test_r_csv_format(tmp_path):
    import pandas as pd
    p = tmp_path / "a.csv"
    rio.write_r_csv(pd.DataFrame({"a": [1.0, float("nan"), 0.1 + 0.2], "b": ["x", "y", "z"]}), str(p))
    assert p.read_text() == '"a","b"\n1,"x"\nNA,"y"\n0.3,"z"\n'
