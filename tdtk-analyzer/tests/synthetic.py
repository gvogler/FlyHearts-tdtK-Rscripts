"""Synthetic test data: a tiny OLE2 (CFB v3) writer and fake tdtK heart movies.

Real .cxd files come from HCImage; this writer produces the same layout
(Field Data/Field n/i_Image1/Bitmap 1 streams + metadata streams) so the
whole pipeline can be exercised without real data.
"""

from __future__ import annotations

import math
import struct

import numpy as np

SECTOR = 512
MINI = 64
CUTOFF = 4096
FREESECT, ENDOFCHAIN, FATSECT, NOSTREAM = 0xFFFFFFFF, 0xFFFFFFFE, 0xFFFFFFFD, 0xFFFFFFFF


def write_cfb(path: str, streams: dict[str, bytes]) -> None:
    """Write an OLE2 compound file. ``streams`` maps 'A/B/name' -> data."""
    # ---- build the storage tree
    root = {"name": "Root Entry", "type": 5, "kids": {}}
    for p, data in streams.items():
        node = root
        parts = p.split("/")
        for part in parts[:-1]:
            node = node["kids"].setdefault(part, {"name": part, "type": 1, "kids": {}})
        node["kids"][parts[-1]] = {"name": parts[-1], "type": 2, "data": data}
    entries = []

    def add(node):
        node["id"] = len(entries)
        entries.append(node)
        kids = sorted(node.get("kids", {}).values(), key=lambda k: (len(k["name"]), k["name"].upper()))
        for k in kids:
            add(k)
        node["child_ids"] = [k["id"] for k in kids]

    add(root)

    # ---- allocate: big streams, mini stream, minifat, directory, fat
    fat: list[int] = []
    chunks: list[bytes] = []

    def alloc(data: bytes) -> int:
        n = max(1, math.ceil(len(data) / SECTOR))
        start = len(fat)
        for i in range(n):
            fat.append(start + i + 1 if i < n - 1 else ENDOFCHAIN)
        chunks.append(data.ljust(n * SECTOR, b"\0"))
        return start

    mini = bytearray()
    minifat: list[int] = []
    for e in entries:
        if e["type"] != 2:
            continue
        d = e["data"]
        if len(d) >= CUTOFF:
            e["start"] = alloc(d)
        else:
            n = max(1, math.ceil(len(d) / MINI))
            e["start"] = len(minifat)
            for i in range(n):
                minifat.append(len(minifat) + 1 if i < n - 1 else ENDOFCHAIN)
            mini += d.ljust(n * MINI, b"\0")
    root["start"] = alloc(bytes(mini)) if mini else ENDOFCHAIN
    root["size"] = len(mini)
    mf_bytes = b"".join(struct.pack("<I", v) for v in minifat)
    mf_start = alloc(mf_bytes) if minifat else ENDOFCHAIN
    n_mf = math.ceil(len(mf_bytes) / SECTOR) if minifat else 0

    def dir_entry(e) -> bytes:
        name = e["name"].encode("utf-16-le") + b"\0\0"
        b = bytearray(128)
        b[0:len(name)] = name
        struct.pack_into("<H", b, 64, len(name))
        b[66] = e["type"]
        b[67] = 1  # black
        struct.pack_into("<III", b, 68, e.get("left", NOSTREAM), e.get("right", NOSTREAM), e.get("tree", NOSTREAM))
        size = len(e["data"]) if e["type"] == 2 else e.get("size", 0)
        start = e.get("start", ENDOFCHAIN if e["type"] != 2 else 0)
        struct.pack_into("<IQ", b, 116, start if e["type"] != 1 else 0, size)
        return bytes(b)

    def balance(ids):  # siblings form a balanced binary search tree (like real files)
        if not ids:
            return NOSTREAM
        mid = len(ids) // 2
        node = entries[ids[mid]]
        node["left"] = balance(ids[:mid])
        node["right"] = balance(ids[mid + 1:])
        return ids[mid]

    for e in entries:
        if e.get("child_ids"):
            e["tree"] = balance(e["child_ids"])
    dirs = b"".join(dir_entry(e) for e in entries)
    dir_start = alloc(dirs)

    n_other = len(fat)
    per = SECTOR // 4
    n_fat, n_difat = 1, 0
    while True:
        n_difat = max(0, math.ceil((n_fat - 109) / (per - 1)))
        if n_fat * per >= n_other + n_fat + n_difat:
            break
        n_fat += 1
    fat_start = len(fat)
    fat.extend([FATSECT] * n_fat)
    difat_start = len(fat)
    fat.extend([0xFFFFFFFC] * n_difat)  # DIFSECT
    fat.extend([FREESECT] * (n_fat * per - len(fat)))
    fat_bytes = b"".join(struct.pack("<I", v) for v in fat)
    fat_ids = [fat_start + i for i in range(n_fat)]
    difat_bytes = b""
    rest = fat_ids[109:]
    for d in range(n_difat):
        part = rest[d * (per - 1):(d + 1) * (per - 1)]
        part = part + [FREESECT] * (per - 1 - len(part))
        nxt = difat_start + d + 1 if d < n_difat - 1 else ENDOFCHAIN
        difat_bytes += struct.pack(f"<{per}I", *part, nxt)

    h = bytearray(SECTOR)
    h[0:8] = bytes.fromhex("D0CF11E0A1B11AE1")
    struct.pack_into("<HHHHH", h, 24, 0x3E, 3, 0xFFFE, 9, 6)
    struct.pack_into("<IIIIIIIII", h, 40, 0, n_fat, dir_start, 0, CUTOFF, mf_start, n_mf,
                     difat_start if n_difat else ENDOFCHAIN, n_difat)
    head = fat_ids[:109] + [FREESECT] * (109 - min(109, n_fat))
    struct.pack_into("<109I", h, 76, *head)
    with open(path, "wb") as fh:
        fh.write(h)
        for c in chunks:
            fh.write(c)
        fh.write(fat_bytes)
        fh.write(difat_bytes)


def heart_movie(T=500, Y=56, X=240, dt=0.005, period=0.25, speed=3000.0, seed=0, reverse=False) -> np.ndarray:
    """uint8 movie [T, Y, X] of a fluorescent heart tube contracting as a travelling wave.

    Bright heart along x in [0, 130) and [160, 235); contractions travel along x
    at ``speed`` px/s (towards smaller x if ``reverse``)."""
    rng = np.random.default_rng(seed)
    t = np.arange(T) * dt
    x = np.arange(X)
    yc = (Y - 1) / 2
    heart = ((x < 130) | ((x >= 160) & (x < 235))).astype(float)
    amp = 0.45 * (1 + 0.08 * np.sin(2 * np.pi * x / 45.0))
    lag = (X - x if reverse else x) / speed
    phase = ((t[:, None] - lag[None, :]) % period) / period          # [T, X]
    contraction = np.where(phase < 0.4, np.sin(np.pi * phase / 0.4) ** 2, 0.0)
    jitter = 1 + 0.15 * rng.standard_normal(T)[:, None]               # beat-to-beat variation
    r_top = 18.0 * (1 - amp[None, :] * contraction * jitter)           # [T, X]
    r_bot = 18.0 * (1 - 0.8 * amp[None, :] * contraction)              # walls move differently
    yy = np.arange(Y)[None, :, None]
    radius = np.where(yy < yc, r_top[:, None, :], r_bot[:, None, :])  # [T, Y, X]
    inside = 1 / (1 + np.exp((np.abs(yy - yc) - radius) / 0.8))        # soft tube
    img = 25 + 190 * inside * heart[None, None, :]
    img += rng.normal(0, 4, img.shape)
    return np.clip(img, 0, 255).astype(np.uint8)


def write_cxd(path: str, movie: np.ndarray, dt: float = 0.005, factor: str = "0.65",
              magnification: str = "1", unix_time: float = 1_669_700_000.0) -> None:
    T, Y, X = movie.shape
    s: dict[str, bytes] = {}
    dbl = lambda v: struct.pack("<d", float(v))  # noqa: E731
    s["File Info/Field Count"] = struct.pack("<i", T)
    s["File Info/File Has Image"] = struct.pack("<h", 1)
    s["File Info/Comments"] = f"factor={factor};um\nmagnification={magnification};\nunits=um;\n".encode()
    s["File Info/Last Field Date & Time"] = dbl(unix_time + 11644444800)
    for i in range(T):
        f = f"Field Data/Field {i + 1}"
        s[f + "/i_Image1/Bitmap 1"] = movie[i].tobytes()
        s[f + "/Details/Time_From_Start"] = dbl(i * dt)
        s[f + "/Details/Time_From_Last"] = dbl(0.0 if i == 0 else dt)
        if i == 0:
            s[f + "/i_Image1/Details/Image_Width"] = dbl(X)
            s[f + "/i_Image1/Details/Image_Height"] = dbl(Y)
            s[f + "/i_Image1/Details/Image_Depth"] = dbl(8 * movie.dtype.itemsize)
    write_cfb(path, s)


def write_imagej_avi(path: str, movie: np.ndarray, fps: float = 7.0) -> None:
    """Uncompressed 8-bit palette AVI as written by ImageJ/Fiji (File > Save As > AVI, no compression):
    bottom-up DIB frames ('00db' chunks), 256-entry grey palette, idx1 index."""
    T, H, W = movie.shape
    stride = (W + 3) // 4 * 4
    frame_bytes = stride * H

    def chunk(fourcc: bytes, data: bytes) -> bytes:
        return fourcc + struct.pack("<I", len(data)) + data + (b"\0" if len(data) % 2 else b"")

    def lst(kind: bytes, data: bytes) -> bytes:
        return b"LIST" + struct.pack("<I", len(data) + 4) + kind + data

    avih = struct.pack("<10I4I", int(round(1e6 / fps)), frame_bytes * int(fps), 0, 0x10, T, 0, 1,
                       frame_bytes, W, H, 0, 0, 0, 0)
    strh = b"vids" + b"DIB " + struct.pack("<IHHIIIIIIIIhhhh", 0, 0, 0, 0, 1000, int(round(fps * 1000)), 0, T,
                                          frame_bytes, 0xFFFFFFFF, 0, 0, 0, W, H)
    bih = struct.pack("<IiiHHIIiiII", 40, W, H, 1, 8, 0, frame_bytes, 0, 0, 256, 0)
    palette = b"".join(bytes((i, i, i, 0)) for i in range(256))
    hdrl = lst(b"hdrl", chunk(b"avih", avih) + lst(b"strl", chunk(b"strh", strh) + chunk(b"strf", bih + palette)))
    frames, index, offset = [], [], 4
    for f in movie:
        rows = np.zeros((H, stride), np.uint8)
        rows[:, :W] = f[::-1]                      # bottom-up
        c = chunk(b"00db", rows.tobytes())
        index.append(b"00db" + struct.pack("<III", 0x10, offset, frame_bytes))
        offset += len(c)
        frames.append(c)
    movi = lst(b"movi", b"".join(frames))
    body = b"AVI " + hdrl + movi + chunk(b"idx1", b"".join(index))
    with open(path, "wb") as fh:
        fh.write(b"RIFF" + struct.pack("<I", len(body)) + body)
