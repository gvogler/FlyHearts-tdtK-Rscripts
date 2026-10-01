import numpy as np

from synthetic import heart_movie, write_cxd
from tdtk_analyzer.cxd import read_cxd_frames, read_cxd_info


def test_cxd_roundtrip(tmp_path):
    movie = heart_movie(T=220, Y=24, X=90)
    p = tmp_path / "A_1_1wf_x.cxd"
    write_cxd(str(p), movie, dt=0.004, factor="1", magnification="2")
    info = read_cxd_info(str(p))
    assert (info.size_x, info.size_y, info.size_t) == (90, 24, 220)
    assert info.pixel_type == "uint8"
    assert info.resolution == 0.65 * 2           # factor 1 -> 0.65 (uncalibrated), times magnification
    assert abs(info.time_interval - 0.004 * 219 / 220) < 1e-12
    np.testing.assert_array_equal(read_cxd_frames(info), movie)
    table = dict(info.meta_table())
    assert list(dict(info.meta_table()))[0] == "sizeX" and len(table) == 21


def test_cxd_uint16(tmp_path):
    movie = (heart_movie(T=210, Y=16, X=64).astype(np.uint16) * 200)
    p = tmp_path / "B_1_1wm_x.cxd"
    write_cxd(str(p), movie)
    np.testing.assert_array_equal(read_cxd_frames(read_cxd_info(str(p))), movie)
