"""Checks whether a movie's metadata is complete and plausible, and parses user corrections.

Every problem is an ``Issue(field, level, message)``:
  level "error": the analysis cannot run (e.g. no frame interval / pixel size)
  level "warn" : the value is suspicious, assumed or inconsistent - please check
  level "info" : worth knowing (e.g. the date comes from the file's modification time)
"""

from __future__ import annotations

import datetime as dt
import re
from dataclasses import dataclass

from .base import RED_CHANNEL_HINTS, ImportOptions, MovieInfo

FIELDS = {"time_interval": "Frame interval", "pixel_size": "Pixel size", "created": "Recording date",
          "channel": "Channel", "frames": "Frames", "axes": "Axes", "file": "File"}
# typical playback rates of exported videos (not the camera frame rate)
PLAYBACK_FPS = (5, 7, 10, 12, 15, 20, 24, 25, 29.97, 30, 50, 59.94, 60)
INTERVAL_RANGE_S = (5e-5, 10.0)       # 20 kHz .. 0.1 Hz
PIXEL_RANGE_UM = (0.02, 30.0)


@dataclass
class Issue:
    field: str
    level: str
    message: str

    def __str__(self):
        return f"{FIELDS.get(self.field, self.field)}: {self.message}"


def check_metadata(m: MovieInfo, opts: ImportOptions) -> list[Issue]:
    out: list[Issue] = []
    src = m.sources

    # ---- frame interval --------------------------------------------------
    if m.time_interval is None or m.time_interval <= 0:
        if opts.default_interval_ms > 0:
            out.append(Issue("time_interval", "info", f"not in the file - using the default {opts.default_interval_ms:g} ms"))
        else:
            out.append(Issue("time_interval", "error", "missing - enter the frame interval (e.g. '5 ms' or '200 fps')"))
    else:
        ti, s = m.time_interval, src.get("time_interval", "file")
        if not INTERVAL_RANGE_S[0] <= ti <= INTERVAL_RANGE_S[1]:
            out.append(Issue("time_interval", "warn", f"{ti * 1000:.4g} ms looks implausible - please check"))
        if s == "video playback rate":
            fps = 1.0 / ti
            if any(abs(fps - p) < 0.02 for p in PLAYBACK_FPS):
                out.append(Issue("time_interval", "warn", f"{fps:.4g} fps is a typical video playback rate, probably "
                                                          "not the camera frame rate - enter the real frame interval"))
            else:
                out.append(Issue("time_interval", "warn", f"taken from the video frame rate ({fps:.4g} fps) - "
                                                          "confirm it is the recording frame rate"))
        elif s.startswith("assumed"):
            out.append(Issue("time_interval", "warn", s))
        tm = m.timing
        if tm:
            if s != "entered" and "time stamps" not in s and abs(tm["median"] - ti) / ti > 0.02:
                out.append(Issue("time_interval", "warn", f"header says {ti * 1000:.4g} ms but the frame time stamps "
                                                          f"give {tm['median'] * 1000:.4g} ms"))
            if tm["cv"] > 0.05:
                out.append(Issue("time_interval", "warn", f"irregular frame timing (variation {tm['cv'] * 100:.0f} %)"))
            if tm["max_gap"] > 1.5:
                out.append(Issue("time_interval", "warn", f"possible dropped frames: largest gap is "
                                                          f"{tm['max_gap']:.1f}x the normal interval (the analysis "
                                                          "assumes evenly spaced frames)"))

    # ---- pixel size ---------------------------------------------------------
    if m.pixel_size is None or m.pixel_size <= 0:
        if opts.default_pixel_um > 0:
            out.append(Issue("pixel_size", "info", f"not in the file - using the default {opts.default_pixel_um:g} µm"))
        else:
            out.append(Issue("pixel_size", "error", "missing - enter the pixel size in µm (camera pixel × binning / "
                                                    "magnification)"))
    else:
        px, s = m.pixel_size, src.get("pixel_size", "file")
        if s.startswith("assumed"):
            out.append(Issue("pixel_size", "warn", s))
        if not PIXEL_RANGE_UM[0] <= px <= PIXEL_RANGE_UM[1]:
            out.append(Issue("pixel_size", "warn", f"{px:.4g} µm per pixel looks implausible - please check"))
        elif px == 1.0 and s != "entered":
            out.append(Issue("pixel_size", "warn", "exactly 1 µm per pixel usually means 'not calibrated'"))

    # ---- recording date ---------------------------------------------------------
    if m.created_unix is None:
        out.append(Issue("created", "warn", "unknown - the tables will have no recording day"))
    else:
        s = src.get("created", "file")
        if s == "file modification time":
            out.append(Issue("created", "info", "taken from the file's modification date - may be a copy date"))
        year = dt.datetime.fromtimestamp(m.created_unix).year
        if year < 1995 or m.created_unix > dt.datetime.now().timestamp() + 86400:
            out.append(Issue("created", "warn", f"{dt.datetime.fromtimestamp(m.created_unix):%Y-%m-%d} looks wrong"))

    # ---- channels / axes ------------------------------------------------------------
    if m.size_c > 1:
        spec = (m.channel if m.channel not in (None, "") else opts.channel) or "auto"
        if spec.lower() == "auto":
            named = [n for n in m.channel_names if any(h in str(n).lower() for h in RED_CHANNEL_HINTS)]
            if not m.channel_names or len(m.channel_names) < m.size_c:
                out.append(Issue("channel", "warn", f"{m.size_c} channels without names - channel 0 is used; "
                                                    "check it in the preview"))
            elif not named:
                out.append(Issue("channel", "warn", f"no tdTomato-like channel among {', '.join(map(str, m.channel_names))} "
                                                    "- channel 0 is used"))
    for n in m.notes:
        fld = "axes" if ("axis" in n or "line scan" in n) else "file"
        level = "warn" if ("used as time" in n or "line scan" in n or "lossy" in n or "differs" in n) else "info"
        out.append(Issue(fld, level, n))
    return out


# --------------------------------------------------------------------------- parsing user input

def parse_interval(text: str) -> float | None:
    """'5', '5 ms', '0.005 s', '200 fps', '200 Hz', '5000 us' -> seconds. Plain numbers are ms."""
    t = text.strip().lower().replace(",", ".").replace("µ", "u")
    if not t:
        return None
    m = re.fullmatch(r"([0-9]*\.?[0-9]+(?:e-?\d+)?)\s*(ms|s|sec|us|min|fps|hz|/s)?", t)
    if not m:
        raise ValueError(f"cannot read a frame interval from '{text}' (examples: 5 ms, 0.005 s, 200 fps)")
    v, unit = float(m.group(1)), m.group(2) or "ms"
    if v <= 0:
        raise ValueError("the frame interval must be positive")
    return {"ms": v / 1000, "s": v, "sec": v, "us": v / 1e6, "min": v * 60}.get(unit, 1.0 / v if unit in ("fps", "hz", "/s") else v)


def parse_pixel(text: str) -> float | None:
    """'0.65', '0.65 um', '650 nm' -> micrometre per pixel. Plain numbers are µm."""
    t = text.strip().lower().replace(",", ".").replace("µ", "u")
    if not t:
        return None
    m = re.fullmatch(r"([0-9]*\.?[0-9]+(?:e-?\d+)?)\s*(um|micron|microns|nm|mm)?", t)
    if not m:
        raise ValueError(f"cannot read a pixel size from '{text}' (examples: 0.65, 0.65 um, 650 nm)")
    v = float(m.group(1)) * {"nm": 1e-3, "mm": 1e3}.get(m.group(2) or "um", 1.0)
    if v <= 0:
        raise ValueError("the pixel size must be positive")
    return v


def parse_date(text: str) -> float | None:
    t = text.strip()
    if not t:
        return None
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d", "%d.%m.%Y %H:%M", "%d.%m.%Y", "%m/%d/%Y"):
        try:
            return dt.datetime.strptime(t, fmt).timestamp()
        except ValueError:
            pass
    raise ValueError(f"cannot read a date from '{text}' (use YYYY-MM-DD or YYYY-MM-DD HH:MM)")


def format_date(ts: float | None) -> str:
    return dt.datetime.fromtimestamp(ts).strftime("%Y-%m-%d %H:%M") if ts else ""
