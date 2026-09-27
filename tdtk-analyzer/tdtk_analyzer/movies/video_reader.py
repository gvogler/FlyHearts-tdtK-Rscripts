"""Ordinary video files (.avi, .mp4, .mov, .mkv) via PyAV (FFmpeg).

Note: most video codecs are lossy and 8-bit; use them only when the microscope
software cannot export a scientific format."""

from __future__ import annotations

import os

import numpy as np

from .base import ImportOptions, MovieInfo


def _gray_format(stream) -> bool:
    return str(stream.codec_context.pix_fmt or "").startswith("gray")


def probe(path: str, rel: str) -> list[MovieInfo]:
    import av

    with av.open(path) as c:
        st = c.streams.video[0]
        fps = float(st.average_rate) if st.average_rate else None
        n = st.frames
        if not n:                              # some containers do not store the frame count
            n = sum(1 for _ in c.decode(st))
        gray = _gray_format(st)
        info = MovieInfo(path=path, rel=rel, format=f"Video ({st.codec_context.name})", reader="video",
                         name=os.path.basename(path), size_t=int(n), size_y=st.height, size_x=st.width,
                         size_c=1 if gray else 3, dtype="uint8", bits=8,
                         time_interval=(1.0 / fps) if fps else None, created_unix=os.path.getmtime(path),
                         channel_names=[] if gray else ["red", "green", "blue"])
    info.notes.append("video file: usually lossy 8-bit compression; pixel size must be entered")
    return [info]


def read(info: MovieInfo, opts: ImportOptions, max_frames: int | None = None) -> np.ndarray:
    import av

    ch = info.channel_index(opts)
    frames = []
    with av.open(info.path) as c:
        st = c.streams.video[0]
        gray = _gray_format(st)
        for i, fr in enumerate(c.decode(st)):
            if max_frames and i >= max_frames:
                break
            if gray:
                frames.append(fr.to_ndarray(format="gray"))
            else:
                frames.append(fr.to_ndarray(format="rgb24")[..., ch])
    return np.stack(frames)
