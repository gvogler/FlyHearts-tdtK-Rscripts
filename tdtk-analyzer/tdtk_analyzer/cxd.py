"""Reader for Hamamatsu HCImage / Compix SimplePCI ``.cxd`` movies.

A .cxd file is an OLE2 compound document. Every frame is stored as a stream
``Field Data/Field <n>/i_Image1/Bitmap 1`` (raw little-endian pixels or a small
TIFF), and the metadata are small streams (8-byte doubles or text). The logic
follows Bio-Formats' ``PCIReader`` (used by the R script via bftools/RBioFormats)
so that sizes, frame order and timing are identical - but no Java is needed.
"""

from __future__ import annotations

import io
import re
import struct
import sys
from dataclasses import dataclass, field

import numpy as np
import olefile

# HCImage stores dates as seconds since 1601-01-01 08:00 - the R script
# subtracts this constant to get Unix time, and so do we.
CXD_EPOCH_OFFSET = 11644444800
MISSING_DATE_VALUE = 11644444801  # used by the R script when the date is missing


class CxdError(Exception):
    pass


@dataclass
class CxdInfo:
    path: str
    size_x: int = 0
    size_y: int = 0
    size_z: int = 1
    size_c: int = 1
    size_t: int = 0
    image_count: int = 0
    bits_per_pixel: int = 0
    pixel_type: str = "uint16"
    dimension_order: str = "XYCTZ"
    time_from_start: dict = field(default_factory=dict)  # field index (1-based) -> s
    time_from_last: dict = field(default_factory=dict)
    last_field_date: float | None = None
    comments: dict = field(default_factory=dict)
    frame_streams: list = field(default_factory=list)  # ordered by plane index
    tiff_frames: bool = False

    # ---- values the R script derives from the metadata -------------------
    @property
    def time_interval(self) -> float:
        vals = np.array(list(self.time_from_last.values()), dtype=float)
        vals = vals[~np.isnan(vals)]
        return float(vals.mean()) if vals.size else float("nan")

    def _comment_number(self, key: str):
        v = self.comments.get(key)
        if v is None:
            return None
        try:
            return float(v.split(";")[0].strip())
        except ValueError:
            return float("nan")

    @property
    def scale_factor(self):
        return self._comment_number("factor")

    @property
    def magnification(self) -> float:
        m = self._comment_number("magnification")
        if m is None or np.isnan(m):
            return 1.0
        return m

    @property
    def resolution(self) -> float:
        """um per pixel, with the R script's corrections."""
        sf = self.scale_factor
        if sf is None or np.isnan(sf):
            raise CxdError(
                f"No calibration data found for file {self.path}. Please move this file into a "
                "different folder and re-run without this file.")
        if sf == 1:
            sf = 0.65  # uncalibrated -> 10x setting for the 6.5 um CMOS chip
        # (the R script's attempt to fix magnification 2 & factor 0.65 had no effect,
        #  so it is intentionally not applied here either)
        return sf * self.magnification

    @property
    def created_unix(self) -> float:
        v = self.last_field_date if self.last_field_date is not None else MISSING_DATE_VALUE
        return v - CXD_EPOCH_OFFSET

    def meta_table(self) -> list[tuple[str, object]]:
        """The 21 values written to '<movie>._new_meta_data.csv' (order matters)."""
        return [
            ("sizeX", self.size_x),
            ("sizeY", self.size_y),
            ("sizeZ", self.size_z),
            ("sizeC", 1),
            ("sizeT", self.size_t),
            ("pixelType", self.pixel_type),
            ("bitsPerPixel", self.bits_per_pixel),
            ("imageCount", self.image_count),
            ("dimensionOrder", self.dimension_order),
            ("orderCertain", self.dimension_order),  # (sic) as in the R script
            ("rgb", "false"),
            ("littleEndian", "true"),
            ("interleaved", "false"),
            ("falseColor", 0),
            ("metadataComplete", 1),
            ("thumbnail", 0),
            ("series", 1),
            ("resolutionLevel", 1),
            ("time_interval", self.time_interval),
            ("resolution", self.resolution),
            ("created_unix_from_file", self.created_unix),
        ]


_FIELD_RE = re.compile(r"Field (\d+)")
_IMAGE_RE = re.compile(r"Image(\d+)")


class _OleFile(olefile.OleFileIO):
    """olefile checks every stream against a *list* of the streams seen so far
    ("stream referenced twice"), which is quadratic in the number of frames; use sets."""

    def _check_duplicate_stream(self, first_sect, minifat=False):
        if not minifat and first_sect in (olefile.DIFSECT, olefile.FATSECT, olefile.ENDOFCHAIN, olefile.FREESECT):
            return
        seen = self.__dict__.setdefault("_seen_streams", (set(), set()))[bool(minifat)]
        if first_sect in seen:
            self._raise_defect(olefile.DEFECT_INCORRECT, "Stream referenced twice")
        seen.add(first_sect)


def _stream_index(ole) -> dict:
    """{path tuple: directory entry} of every stream, in the order of ``ole.listdir()``.

    olefile finds a stream by searching the directory for its path on every
    ``openstream``/``get_size`` call, so reading all frames of a movie took time
    proportional to frames² (about 10 s just for the metadata of 6000 frames).
    One index per file makes every lookup a dictionary access."""
    out = {}

    def walk(node, prefix):
        for kid in node.kids:
            if kid.entry_type == olefile.STGTY_STORAGE:
                walk(kid, prefix + (kid.name,))
            elif kid.entry_type == olefile.STGTY_STREAM:
                out[prefix + (kid.name,)] = kid

    walk(ole.root, ())
    return out


def _read(ole, de, n: int = -1) -> bytes:
    """Contents (or the first ``n`` bytes) of the stream with directory entry ``de``."""
    return ole._open(de.isectStart, de.size).read(n)


def _read_double(ole, de) -> float:
    return struct.unpack("<d", _read(ole, de, 8))[0]


def read_cxd_info(path: str) -> CxdInfo:
    """Parse the metadata and the frame layout of a .cxd file."""
    if not olefile.isOleFile(path):
        raise CxdError(f"{path} is not a CXD (OLE2) file")
    # olefile walks the directory tree recursively; movies have thousands of frames
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 20000))
    info = CxdInfo(path=path)
    frames: list[tuple[str, list[str]]] = []
    unique_z: list[float] = []
    first_z = second_z = 0.0
    mode = 0
    group_selected = None
    with _OleFile(path) as ole:
        index = _stream_index(ole)
        if not index:
            raise CxdError("No files were found - the .cxd may be corrupt.")
        for entry, de in index.items():
            rel = entry[-1].strip()
            parent = "/".join(entry[:-1])
            is_bitmap = rel.startswith("Bitmap") or (rel == "Data" and "Image" in parent)
            if is_bitmap:
                frames.append((parent, list(entry)))
                continue
            size = de.size
            if rel == "Field Count":
                info.image_count = struct.unpack("<i", _read(ole, de, 4))[0]
            elif rel == "File Has Image":
                if struct.unpack("<h", _read(ole, de, 2))[0] == 0:
                    raise CxdError("This file does not contain image data.")
            elif "Image_Depth" in rel:
                bits = int(_read_double(ole, de))
                info.bits_per_pixel = bits
                while bits % 8 != 0 or bits == 0:
                    bits += 1
                if bits % 3 == 0:
                    info.size_c = 3
                    bits //= 3
                    info.bits_per_pixel //= 3
                nbytes = bits // 8
                info.pixel_type = {1: "uint8", 2: "uint16", 4: "uint32"}.get(nbytes, "uint16")
            elif "Image_Height" in rel and info.size_y == 0:
                info.size_y = int(_read_double(ole, de))
            elif "Image_Width" in rel and info.size_x == 0:
                info.size_x = int(_read_double(ole, de))
            elif "Time_From_Start" in rel or "Time_From_Last" in rel:
                m = _FIELD_RE.findall(parent)
                if m and size >= 8:
                    target = info.time_from_start if "Time_From_Start" in rel else info.time_from_last
                    target[int(m[-1])] = _read_double(ole, de)
            elif rel.endswith("Position_Z") and size >= 8:
                z = _read_double(ole, de)
                if z not in unique_z:
                    unique_z.append(z)
                if "Field 1/" in parent + "/":
                    first_z = z
                elif "Field 2/" in parent + "/":
                    second_z = z
            elif "Last Field Date" in rel and size >= 8:
                info.last_field_date = _read_double(ole, de)
            elif rel == "GroupMode":
                mode = struct.unpack("<i", _read(ole, de, 4))[0]
            elif rel == "GroupSelectedFields":
                group_selected = size // 8
            elif rel == "Comments":
                text = _read(ole, de).decode("latin-1", errors="replace")
                for line in text.replace("\r", "\n").split("\n"):
                    if "=" in line:
                        k, v = line.split("=", 1)
                        info.comments[k.strip()] = v.strip()

        # ---- Z / T bookkeeping (as PCIReader) ----------------------------
        size_z = group_selected if (group_selected and mode != 0) else 0
        if size_z <= 1 or (info.image_count % size_z) != 0:
            size_z = len(unique_z) if unique_z else 1
        if info.image_count == 0:
            info.image_count = len(frames)
        size_t = info.image_count // size_z
        while size_z * size_t < info.image_count:
            size_z += 1
            size_t = info.image_count // size_z
        info.size_z, info.size_t = size_z, size_t
        info.dimension_order = "XYCZT" if abs(first_z - second_z) > 1e-6 else "XYCTZ"
        info.image_count = size_z * size_t

        # ---- order frames by field number --------------------------------
        def plane_index(parent: str) -> int:
            f = _FIELD_RE.findall(parent)
            im = _IMAGE_RE.findall(parent)
            field_no = int(f[-1]) if f else 1
            image_no = int(im[-1]) if im else 1
            return info.size_c * (field_no - 1) + (image_no - 1)

        frames.sort(key=lambda fr: plane_index(fr[0]))
        info.frame_streams = [e for _, e in frames]
        if not info.frame_streams:
            raise CxdError("This file does not contain image data.")

        # ---- padded rows: widen X like Bio-Formats does -------------------
        bpp = np.dtype(info.pixel_type).itemsize
        first_de = index[tuple(info.frame_streams[0])]
        first = _read(ole, first_de, 16)
        info.tiff_frames = first[:4] in (b"II*\x00", b"MM\x00*")
        if not info.tiff_frames:
            length = first_de.size
            expected = info.size_x * info.size_y * bpp * info.size_c
            if length > expected and info.size_y > 0:
                extra = (length - expected) // (info.size_y * bpp * info.size_c)
                info.size_x += int(extra)
    return info


def read_cxd_frames(info: CxdInfo, n_frames: int | None = None) -> np.ndarray:
    """Read the movie as an array [T, Y, X] (native integer type).

    Frames beyond the planes listed in the metadata are ignored, as in
    Bio-Formats. ``n_frames`` limits the number of frames read."""
    import tifffile

    n = min(len(info.frame_streams), info.image_count or len(info.frame_streams))
    if n_frames is not None:
        n = min(n, n_frames)
    dtype = np.dtype(info.pixel_type).newbyteorder("<")
    plane = info.size_x * info.size_y
    out = np.empty((n, info.size_y, info.size_x), dtype=dtype.newbyteorder("="))
    with _OleFile(info.path) as ole:
        index = _stream_index(ole)
        for i in range(n):
            data = _read(ole, index[tuple(info.frame_streams[i])])
            if info.tiff_frames:
                out[i] = tifffile.imread(io.BytesIO(data))
            else:
                out[i] = np.frombuffer(data, dtype=dtype, count=plane).reshape(info.size_y, info.size_x)
    return out
