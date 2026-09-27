"""Import of different microscope formats (files written with the formats' own libraries)."""

import numpy as np
import pytest
import tifffile

from synthetic import heart_movie, write_cxd
from tdtk_analyzer.movies import (ImportOptions, apply_table, find_movie_files, load_movie, movie_status,
                                  probe_file, scan_movies, to_table)
from tdtk_analyzer.movies.base import auto_rotation

MOVIE = heart_movie(T=220, Y=32, X=120, seed=3)          # uint8 [T, Y, X]


def test_ome_tiff_two_channels_picks_tdtomato(tmp_path):
    p = tmp_path / "A_1_1wf.ome.tif"
    green = (255 - MOVIE) // 4
    data = np.stack([green, MOVIE], axis=1).astype(np.uint16)       # T C Y X
    tifffile.imwrite(p, data, ome=True, metadata={
        "axes": "TCYX", "PhysicalSizeX": 0.8, "PhysicalSizeXUnit": "µm", "TimeIncrement": 5, "TimeIncrementUnit": "ms",
        "Channel": {"Name": ["GFP", "tdTomato"]}})
    (m,) = probe_file(str(p))
    assert m.format == "OME-TIFF" and (m.size_t, m.size_c, m.size_y, m.size_x) == (220, 2, 32, 120)
    assert m.pixel_size == pytest.approx(0.8) and m.time_interval == pytest.approx(0.005)
    assert m.channel_names == ["GFP", "tdTomato"]
    movie, rot = load_movie(m, ImportOptions())
    np.testing.assert_array_equal(movie, MOVIE)                     # 'auto' chose tdTomato
    g, _ = load_movie(m, ImportOptions(channel="GFP"))
    np.testing.assert_array_equal(g, green)
    assert load_movie(m, ImportOptions(), max_frames=10)[0].shape == (10, 32, 120)


def test_imagej_hyperstack(tmp_path):
    p = tmp_path / "B_2_1wm.tif"
    tifffile.imwrite(p, MOVIE, imagej=True, resolution=(1 / 0.65, 1 / 0.65),
                     metadata={"axes": "TYX", "finterval": 0.004, "unit": "um"})
    (m,) = probe_file(str(p))
    assert m.format == "ImageJ TIFF"
    assert m.time_interval == pytest.approx(0.004) and m.pixel_size == pytest.approx(0.65)
    np.testing.assert_array_equal(load_movie(m, ImportOptions())[0], MOVIE)


def test_stack_saved_as_z_uses_z_as_time(tmp_path):
    p = tmp_path / "C_1_1wf.tif"
    tifffile.imwrite(p, MOVIE, imagej=True, metadata={"axes": "ZYX"})
    (m,) = probe_file(str(p))
    assert m.size_t == 220 and any("used as time" in n for n in m.notes)
    np.testing.assert_array_equal(load_movie(m, ImportOptions())[0], MOVIE)


def test_plain_tiff_needs_interval_and_pixel(tmp_path):
    p = tmp_path / "D_1_1wf.tiff"
    tifffile.imwrite(p, MOVIE)
    (m,) = probe_file(str(p))
    assert m.time_interval is None
    level, msg = movie_status(m, ImportOptions())
    assert level == "error" and "frame interval" in msg
    opts = ImportOptions(default_interval_ms=5, default_pixel_um=0.65)
    assert movie_status(m, opts)[0] in ("ok", "warn")      # warn: 'image axis used as time'
    # or per movie through the import table
    t = to_table([m])
    t.loc[0, "frame_interval_ms"], t.loc[0, "pixel_size_um"], t.loc[0, "name"] = 4, 1.3, "D_9_1wf"
    apply_table([m], t.astype(str))
    assert m.time_interval == pytest.approx(0.004) and m.pixel_size == pytest.approx(1.3)
    assert m.name == "D_9_1wf.tiff"                        # the extension is kept


def test_czi(tmp_path):
    from pylibCZIrw import czi as pyczi

    p = tmp_path / "E_1_1wf.czi"
    with pyczi.create_czi(str(p), exist_ok=True) as w:
        for t in range(MOVIE.shape[0]):
            w.write(data=MOVIE[t].astype(np.uint16)[..., None], plane={"T": t, "C": 0, "Z": 0})
        w.write_metadata(document_name="test", channel_names={0: "tdTomato"}, scale_x=0.8e-6, scale_y=0.8e-6)
    (m,) = probe_file(str(p))
    assert m.format == "Zeiss CZI" and (m.size_t, m.size_y, m.size_x) == (220, 32, 120)
    assert m.pixel_size == pytest.approx(0.8)
    np.testing.assert_array_equal(load_movie(m, ImportOptions())[0], MOVIE)


def test_video(tmp_path):
    import av

    p = tmp_path / "F_1_1wf.mp4"
    frames = MOVIE[:, :, :112]                                   # even width for the codec
    with av.open(str(p), "w") as out:
        st = out.add_stream("libx264", rate=200)
        st.width, st.height, st.pix_fmt = 112, 32, "yuv420p"
        st.options = {"crf": "10"}
        for f in frames:
            for pkt in st.encode(av.VideoFrame.from_ndarray(np.repeat(f[..., None], 3, axis=2), format="rgb24")):
                out.mux(pkt)
        for pkt in st.encode():
            out.mux(pkt)
    (m,) = probe_file(str(p))
    assert m.size_t == 220 and m.time_interval == pytest.approx(1 / 200)
    movie, _ = load_movie(m, ImportOptions())
    assert movie.shape == (220, 32, 112)
    assert np.corrcoef(movie.ravel().astype(float), frames.ravel().astype(float))[0, 1] > 0.95


def test_auto_rotation_detects_vertical_heart():
    assert auto_rotation(MOVIE) == "0"
    assert auto_rotation(np.transpose(MOVIE, (0, 2, 1))) == "90"


def test_scan_skips_generated_files(tmp_path):
    write_cxd(str(tmp_path / "G_1_1wf.cxd"), MOVIE)
    tifffile.imwrite(tmp_path / "G_1_1wf.cxd_peak_1_at Xpos_10.tiff", MOVIE[0])
    tifffile.imwrite(tmp_path / "G_1_1wf._SD and peaklines.tiff", MOVIE[0])
    (tmp_path / "notes.txt").write_text("x")
    assert find_movie_files(str(tmp_path)) == ["G_1_1wf.cxd"]
    (m,) = scan_movies(str(tmp_path))
    assert m.format == "Hamamatsu CXD" and m.name == "G_1_1wf.cxd"


def test_table_keeps_use_choice_after_a_failed_scan(tmp_path):
    from tdtk_analyzer.movies import MovieInfo, save_table, load_table

    bad = MovieInfo(path=str(tmp_path / "x.tif"), rel="x.tif", format="TIFF", reader="tiff", name="x.tif",
                    error="could not read")
    save_table([bad], str(tmp_path))
    good = MovieInfo(path=str(tmp_path / "x.tif"), rel="x.tif", format="TIFF", reader="tiff", name="x.tif")
    apply_table([good], load_table(str(tmp_path)))
    assert good.use
