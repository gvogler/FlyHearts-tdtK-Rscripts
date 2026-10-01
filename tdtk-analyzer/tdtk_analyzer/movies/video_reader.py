"""Video files: AVI (uncompressed, ImageJ/Fiji 8-bit palette, MJPEG, PNG, FFV1 incl. 16-bit, ...),
MP4, MOV, MKV - via PyAV (FFmpeg).

Videos store a *playback* frame rate, which is often not the camera frame rate
(e.g. Fiji's AVI export defaults to 7 fps), and never a pixel size. Both are
therefore flagged for checking in the import step.
"""

from __future__ import annotations

import os

import numpy as np

from .base import ImportOptions, MovieInfo

GRAY16 = ("gray16le", "gray16be", "gray12le", "gray10le", "gray14le")


def _kind(pix_fmt: str) -> str:
    pix_fmt = pix_fmt or ""
    if pix_fmt in GRAY16 or (pix_fmt.startswith("gray") and ("16" in pix_fmt or "12" in pix_fmt or "10" in pix_fmt)):
        return "gray16"
    if pix_fmt.startswith("gray") or pix_fmt.startswith("ya8"):
        return "gray8"
    return "color"                      # rgb/bgr/yuv/pal8: decoded to RGB


def _is_yuv(pix_fmt: str) -> bool:
    return str(pix_fmt or "").startswith(("yuv", "yuvj", "nv", "uyvy", "yuyv"))


def _frame(fr, kind: str) -> np.ndarray:
    if kind == "gray16":
        return fr.to_ndarray(format="gray16le")
    if kind == "gray8":
        return fr.to_ndarray(format="gray")
    return fr.to_ndarray(format="rgb24")


def probe(path: str, rel: str) -> list[MovieInfo]:
    import av

    with av.open(path) as c:
        st = c.streams.video[0]
        codec = st.codec_context.name
        kind = _kind(st.codec_context.pix_fmt)
        rate = st.average_rate or st.guessed_rate or st.base_rate
        fps = float(rate) if rate else None
        tb = float(st.time_base) if st.time_base else None
        if c.format.name == "avi" and st.frames:
            # AVI time stamps are only frame numbers (no real timing) and collecting them reads
            # the whole file - the frame count in the header is all that is needed
            pts = []
        else:
            # time stamps and frame count from the container packets (no decoding)
            pts = [p.pts for p in c.demux(st) if p.pts is not None and p.size > 0]
            c.seek(0)
        n = st.frames or len(pts)
        fr0 = next(c.decode(st))
        first = _frame(fr0, kind)
        yuv = _is_yuv(st.codec_context.pix_fmt)
        gray_as_color = False
        if kind == "color" and yuv:             # grey content: colour planes stay neutral (128)
            planes = fr0.to_ndarray(format="yuv444p").astype(np.int16)[1:]
            flat = planes.reshape(2, -1)
            centre = np.median(flat, axis=1)
            chroma = np.abs(flat - centre[:, None])        # constant chroma = no colour information
            gray_as_color = bool(np.all(np.abs(centre - 128) <= 3) and np.percentile(chroma, 99.9) <= 4
                                 and chroma.mean() < 1.0)
        elif kind == "color":                   # RGB / palette: all channels equal
            gray_as_color = bool(np.array_equal(first[..., 0], first[..., 1])
                                 and np.array_equal(first[..., 1], first[..., 2]))
    size_c = 3 if (kind == "color" and not gray_as_color) else 1
    info = MovieInfo(path=path, rel=rel, format=f"Video {os.path.splitext(path)[1].lstrip('.').upper()} ({codec})",
                     reader="video", name=os.path.basename(path), size_t=int(n), size_y=first.shape[0],
                     size_x=first.shape[1], size_c=size_c, dtype="uint16" if kind == "gray16" else "uint8",
                     bits=16 if kind == "gray16" else 8, time_interval=(1.0 / fps) if fps else None,
                     created_unix=os.path.getmtime(path),
                     channel_names=["red", "green", "blue"] if size_c == 3 else [],
                     extra={"kind": kind, "gray_as_color": gray_as_color, "yuv": yuv})
    info.sources = {"created": "file modification time"}
    if fps:
        info.sources["time_interval"] = "video playback rate"
    if tb and len(pts) > 2:
        info.set_timing(sorted(p * tb for p in pts))
    if codec not in ("rawvideo", "ffv1", "png", "huffyuv", "utvideo", "ffvhuff", "qtrle", "r210", "v210"):
        info.notes.append(f"lossy video compression ({codec}) - prefer uncompressed or TIFF exports")
    if gray_as_color:
        info.notes.append("grayscale movie stored as color - read as one channel")
    return [info]


def read(info: MovieInfo, opts: ImportOptions, max_frames: int | None = None) -> np.ndarray:
    import av

    kind = info.extra.get("kind", "color")
    ch = info.channel_index(opts)
    if info.extra.get("gray_as_color"):
        if info.extra.get("yuv"):
            kind = "gray8"                     # the luma plane is the grey image
        else:
            ch = 0                             # RGB/palette: exact value of any channel
    frames = []
    with av.open(info.path) as c:
        st = c.streams.video[0]
        for i, fr in enumerate(c.decode(st)):
            if max_frames and i >= max_frames:
                break
            a = _frame(fr, kind)
            frames.append(a[..., ch] if a.ndim == 3 else a)
    return np.stack(frames)
