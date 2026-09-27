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


# ----------------------------------------------------------------------------- AVI

def _write_avi(path, frames, codec, pix, fmt_in, rate=200):
    import av
    with av.open(str(path), "w") as out:
        st = out.add_stream(codec, rate=rate)
        st.width, st.height, st.pix_fmt = frames.shape[2], frames.shape[1], pix
        for f in frames:
            vf = av.VideoFrame.from_ndarray(f if fmt_in != "rgb24" else np.repeat(f[..., None], 3, 2), format=fmt_in)
            for pkt in st.encode(vf if fmt_in == pix else vf.reformat(format=pix)):
                out.mux(pkt)
        for pkt in st.encode():
            out.mux(pkt)


@pytest.mark.parametrize("codec,pix,fmt_in", [("rawvideo", "gray", "gray"), ("rawvideo", "bgr24", "rgb24"),
                                             ("mjpeg", "yuvj420p", "rgb24"), ("png", "rgb24", "rgb24")])
def test_avi_variants(tmp_path, codec, pix, fmt_in):
    frames = MOVIE[:, :, :112]
    p = tmp_path / f"H_1_1wf_{codec}.avi"
    _write_avi(p, frames, codec, pix, fmt_in)
    (m,) = probe_file(str(p))
    assert m.format.startswith("Video AVI") and (m.size_t, m.size_y, m.size_x) == (220, 32, 112)
    assert m.size_c == 1                                      # grey video, even when stored as color
    movie, _ = load_movie(m, ImportOptions())
    if codec in ("rawvideo", "png"):                        # lossless: identical
        np.testing.assert_array_equal(movie, frames)
    else:                                                   # lossy JPEG: close
        assert np.corrcoef(movie.ravel().astype(float), frames.ravel().astype(float))[0, 1] > 0.98


def test_avi_16bit_ffv1(tmp_path):
    frames = MOVIE[:, :, :112].astype(np.uint16) * 257
    p = tmp_path / "I_1_1wf.avi"
    _write_avi(p, frames, "ffv1", "gray16le", "gray16le")
    (m,) = probe_file(str(p))
    assert m.bits == 16
    movie, _ = load_movie(m, ImportOptions())
    assert movie.dtype == np.uint16
    np.testing.assert_array_equal(movie, frames)


def test_imagej_palette_avi_and_playback_rate_warning(tmp_path):
    from synthetic import write_imagej_avi
    from tdtk_analyzer.movies import movie_issues

    p = tmp_path / "J_1_1wf.avi"
    write_imagej_avi(str(p), MOVIE, fps=7)                 # Fiji's default export rate
    (m,) = probe_file(str(p))
    assert (m.size_t, m.size_y, m.size_x, m.size_c) == (220, 32, 120, 1)
    np.testing.assert_array_equal(load_movie(m, ImportOptions())[0], MOVIE)
    msgs = " | ".join(str(i) for i in movie_issues(m, ImportOptions()))
    assert "playback rate" in msgs and "Pixel size: missing" in msgs
    level, _ = movie_status(m, ImportOptions())
    assert level == "error"                               # no pixel size
    m.set_value("time_interval", 0.005)
    m.set_value("pixel_size", 0.65)
    assert movie_status(m, ImportOptions())[0] in ("ok", "warn")
    assert "playback" not in " ".join(str(i) for i in movie_issues(m, ImportOptions()))


# ----------------------------------------------------------------------------- metadata checks

def test_metadata_checks_and_parsing():
    from tdtk_analyzer.movies import MovieInfo, movie_issues
    from tdtk_analyzer.movies.metadata import parse_date, parse_interval, parse_pixel

    assert parse_interval("5") == pytest.approx(0.005)
    assert parse_interval("200 fps") == pytest.approx(0.005)
    assert parse_interval("0,004 s") == pytest.approx(0.004)
    assert parse_interval("5000 µs") == pytest.approx(0.005)
    assert parse_pixel("650 nm") == pytest.approx(0.65)
    assert parse_pixel("0.65 um") == pytest.approx(0.65)
    with pytest.raises(ValueError):
        parse_interval("fast")
    assert parse_date("2024-03-01 10:15") > 0

    m = MovieInfo(path="x.tif", rel="x.tif", format="TIFF", reader="tiff", name="x.tif", size_t=500, size_c=2,
                  time_interval=0.005, pixel_size=1.0, created_unix=None,
                  sources={"time_interval": "file", "pixel_size": "file"})
    m.set_timing(list(np.arange(500) * 0.005) [:250] + list(250 * 0.005 + 0.02 + np.arange(250) * 0.005))
    text = " | ".join(str(i) for i in movie_issues(m, ImportOptions()))
    assert "dropped frames" in text                       # a 4x gap in the time stamps
    assert "not calibrated" in text                       # exactly 1 µm/pixel
    assert "Recording date: unknown" in text
    assert "without names" in text                        # 2 unnamed channels
    m.set_value("pixel_size", 30000.0)
    assert "implausible" in " ".join(str(i) for i in movie_issues(m, ImportOptions()))


def test_cxd_sources_and_timing(tmp_path):
    p = tmp_path / "K_1_1wf.cxd"
    write_cxd(str(p), MOVIE, factor="1")
    (m,) = probe_file(str(p))
    assert m.sources["time_interval"] == "file time stamps" and m.timing["cv"] < 0.01
    assert m.sources["pixel_size"].startswith("assumed")  # factor = 1 -> 0.65 um as in R
    from tdtk_analyzer.movies import movie_issues
    assert any(i.field == "pixel_size" and i.level == "warn" for i in movie_issues(m, ImportOptions()))


def test_user_values_survive_rescan(tmp_path):
    from tdtk_analyzer.movies import load_table, save_table

    p = tmp_path / "L_1_1wf.tif"
    tifffile.imwrite(p, MOVIE)
    (m,) = probe_file(str(p))
    m.set_value("time_interval", 0.004)
    m.set_value("pixel_size", 0.8)
    m.set_value("created_unix", 1_700_000_000.0)
    save_table([m], str(tmp_path))
    (fresh,) = probe_file(str(p))
    apply_table([fresh], load_table(str(tmp_path)))
    assert fresh.time_interval == pytest.approx(0.004) and fresh.pixel_size == pytest.approx(0.8)
    assert fresh.sources["time_interval"] == "entered" and fresh.created_unix is not None
    fresh.reset_to_file()
    assert fresh.time_interval is None and "time_interval" not in fresh.sources
