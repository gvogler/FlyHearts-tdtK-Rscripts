"""TIFF family: OME-TIFF, ImageJ hyperstacks, Micro-Manager, Zeiss LSM, MetaMorph STK, plain stacks."""

from __future__ import annotations

import datetime as dt
import os
import xml.etree.ElementTree as ET

import numpy as np

from .base import MovieInfo, ImportOptions, output_name, pick_time_axis, select_planes

LENGTH_TO_UM = {"m": 1e6, "cm": 1e4, "mm": 1e3, "µm": 1.0, "um": 1.0, "micron": 1.0, "microns": 1.0,
                "micrometer": 1.0, "nm": 1e-3, "å": 1e-4, "a": 1e-4, "inch": 25400.0, "in": 25400.0}
TIME_TO_S = {"s": 1.0, "sec": 1.0, "ms": 1e-3, "µs": 1e-6, "us": 1e-6, "min": 60.0, "h": 3600.0, "hr": 3600.0}


def to_um(value, unit="µm"):
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    f = LENGTH_TO_UM.get(str(unit or "µm").strip().lower().replace("μ", "µ"))
    return v * f if f and v > 0 else None


def to_s(value, unit="s"):
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    f = TIME_TO_S.get(str(unit or "s").strip().lower())
    return v * f if f and v > 0 else None


def parse_iso(s) -> float | None:
    if not s:
        return None
    try:
        return dt.datetime.fromisoformat(str(s).replace("Z", "+00:00")).timestamp()
    except ValueError:
        return None


# --------------------------------------------------------------------------- OME-XML

def parse_ome(xml: str) -> list[dict]:  # noqa: C901
    """Per image: name, sizes, dimension order, pixel size (um), frame interval (s), channels, date."""
    root = ET.fromstring(xml)
    ns = {"ome": root.tag.split("}")[0].strip("{")} if root.tag.startswith("{") else {}
    q = (lambda t: f"ome:{t}") if ns else (lambda t: t)
    images = []
    for img in root.findall(q("Image"), ns):
        px = img.find(q("Pixels"), ns)
        if px is None:
            continue
        a = px.attrib
        sizes = {k: int(a.get(f"Size{k}", 1)) for k in "XYZCT"}
        chans = [c.attrib.get("Name") or c.attrib.get("Fluor") or "" for c in px.findall(q("Channel"), ns)]
        if not any(chans):
            chans = []                        # only IDs ('Channel:0:0') - no real names
        interval = to_s(a.get("TimeIncrement"), a.get("TimeIncrementUnit", "s"))
        deltas = {}
        for pl in px.findall(q("Plane"), ns):
            p = pl.attrib
            if int(p.get("TheC", 0)) == 0 and int(p.get("TheZ", 0)) == 0 and p.get("DeltaT") is not None:
                try:
                    deltas[int(p.get("TheT", 0))] = float(p["DeltaT"]) * TIME_TO_S.get(p.get("DeltaTUnit", "s"), 1.0)
                except ValueError:
                    pass
        stamps = [deltas[k] for k in sorted(deltas)] if len(deltas) > 1 else []
        header_interval = interval
        if stamps and stamps[-1] > stamps[0]:
            d = np.diff(stamps)
            d = d[d > 0]
            # typical spacing: robust against dropped frames (reported separately)
            interval = float(np.median(d)) if d.size else (stamps[-1] - stamps[0]) / (len(stamps) - 1)
        acq = img.find(q("AcquisitionDate"), ns)
        images.append({
            "name": img.attrib.get("Name", ""), "sizes": sizes, "order": a.get("DimensionOrder", "XYCZT"),
            "pixel_size": to_um(a.get("PhysicalSizeX"), a.get("PhysicalSizeXUnit", "µm")),
            "time_interval": interval, "channels": chans, "type": a.get("Type", ""),
            "stamps": stamps, "header_interval": header_interval,
            "created": parse_iso(acq.text if acq is not None else None),
        })
    return images


# --------------------------------------------------------------------------- probe

def _resolution_um(page) -> float | None:
    try:
        xres = page.tags["XResolution"].value
        unit = page.tags["ResolutionUnit"].value if "ResolutionUnit" in page.tags else 2
        ppu = xres[0] / xres[1] if isinstance(xres, tuple) else float(xres)
        if ppu <= 0 or ppu == 1:
            return None
        if int(unit) == 2 and round(ppu) in (72, 96, 150, 300, 600):   # screen/print defaults, not a calibration
            return None
        per = {2: 25400.0, 3: 10000.0}.get(int(unit))
        return per / ppu if per else None
    except Exception:
        return None


def probe(path: str, rel: str) -> list[MovieInfo]:
    import tifffile

    out = []
    fname = os.path.basename(path)
    with tifffile.TiffFile(path) as tf:
        ome = parse_ome(tf.ome_metadata) if tf.is_ome and tf.ome_metadata else []
        fmt = ("OME-TIFF" if tf.is_ome else "Zeiss LSM" if tf.is_lsm else "Micro-Manager TIFF" if tf.is_micromanager
               else "ImageJ TIFF" if tf.is_imagej else "MetaMorph STK" if tf.is_stk else "TIFF stack")
        series = list(tf.series)
        for si, s in enumerate(series):
            axes, shape = s.axes.upper(), tuple(s.shape)
            if "Y" not in axes or "X" not in axes:
                continue
            t_ax, note = pick_time_axis(axes, shape)
            size = lambda a: shape[axes.index(a)] if a in axes else 1  # noqa: E731
            cax = "C" if "C" in axes else ("S" if "S" in axes else "")
            info = MovieInfo(path=path, rel=rel, format=fmt, reader="tiff", series=si,
                             n_series=len(series), size_t=size(t_ax), size_y=size("Y"), size_x=size("X"),
                             size_c=size(cax) if cax else 1, size_z=size("Z") if t_ax != "Z" else 1,
                             dtype=str(s.dtype), bits=np.dtype(s.dtype).itemsize * 8,
                             extra={"axes": axes, "t_axis": t_ax})
            if note:
                info.notes.append(note)
            page = s.pages[0] if len(s.pages) else tf.pages[0]
            info.pixel_size = _resolution_um(page)
            if info.pixel_size:
                info.sources["pixel_size"] = "file (TIFF resolution tag)"
            if cax == "S" and info.size_c == 3:
                info.channel_names = ["red", "green", "blue"]

            if si < len(ome):
                o = ome[si]
                info.series_name = o["name"]
                if o["pixel_size"]:
                    info.pixel_size, info.sources["pixel_size"] = o["pixel_size"], "file (OME)"
                if o["time_interval"]:
                    info.time_interval = o["time_interval"]
                    info.sources["time_interval"] = "file time stamps (OME)" if o["stamps"] else "file (OME)"
                if o["stamps"]:
                    info.set_timing(o["stamps"])
                    if o["header_interval"] and abs(o["header_interval"] - info.time_interval) / info.time_interval > 0.02:
                        info.notes.append(f"OME TimeIncrement {o['header_interval'] * 1000:.4g} ms differs from the "
                                          f"plane time stamps ({info.time_interval * 1000:.4g} ms)")
                info.channel_names = o["channels"] or info.channel_names
                if o["created"]:
                    info.created_unix, info.sources["created"] = o["created"], "file (OME)"
            if tf.is_imagej and tf.imagej_metadata:
                m = tf.imagej_metadata
                if m.get("finterval"):
                    info.time_interval = float(m["finterval"])
                    info.sources["time_interval"] = "file (ImageJ frame interval)"
                elif m.get("fps"):
                    info.time_interval = 1.0 / float(m["fps"])
                    info.sources["time_interval"] = "video playback rate"
                unit = str(m.get("unit", "")).replace("\\u00B5", "µ")
                if unit and unit.lower() not in ("pixel", "pixels", ""):
                    ppu = page.tags["XResolution"].value if "XResolution" in page.tags else None
                    if ppu and to_um(ppu[1] / ppu[0], unit):
                        info.pixel_size = to_um(ppu[1] / ppu[0], unit)
                        info.sources["pixel_size"] = "file (ImageJ calibration)"
                if m.get("Labels") and not info.channel_names and info.size_c > 1:
                    info.channel_names = list(m["Labels"])[: info.size_c]
            if tf.is_micromanager and tf.micromanager_metadata:
                summ = tf.micromanager_metadata.get("Summary", {}) or {}
                if summ.get("Interval_ms"):
                    info.time_interval = float(summ["Interval_ms"]) / 1000.0
                    info.sources["time_interval"] = "file (Micro-Manager interval setting)"
                if summ.get("PixelSize_um"):
                    info.pixel_size = float(summ["PixelSize_um"])
                    info.sources["pixel_size"] = "file (Micro-Manager)"
                if summ.get("ChNames"):
                    info.channel_names = list(summ["ChNames"])
            if tf.is_lsm and tf.lsm_metadata:
                m = tf.lsm_metadata
                if m.get("TimeIntervall"):
                    info.time_interval = float(m["TimeIntervall"])
                    info.sources["time_interval"] = "file (LSM)"
                ts = m.get("TimeStamps")
                if ts is not None and len(ts) > 1:
                    info.time_interval = float((ts[-1] - ts[0]) / (len(ts) - 1))
                    info.sources["time_interval"] = "file time stamps (LSM)"
                    info.set_timing(list(ts))
                if m.get("VoxelSizeX"):
                    info.pixel_size = float(m["VoxelSizeX"]) * 1e6
                    info.sources["pixel_size"] = "file (LSM)"
            if tf.is_stk and tf.stk_metadata:
                m = tf.stk_metadata
                if m.get("XCalibration"):
                    info.pixel_size = float(m["XCalibration"])
                    info.sources["pixel_size"] = "file (MetaMorph)"
                tc = m.get("TimeCreated")
                if tc is not None and len(tc) > 1:
                    try:
                        secs = np.array([t.timestamp() for t in tc])
                        info.time_interval = float((secs[-1] - secs[0]) / (len(secs) - 1)) or None
                        info.sources["time_interval"] = "file time stamps (MetaMorph)"
                        info.set_timing(list(secs))
                    except Exception:
                        pass
            if info.created_unix is None:
                dtag = page.tags.get("DateTime")
                if dtag:
                    try:
                        info.created_unix = dt.datetime.strptime(str(dtag.value), "%Y:%m:%d %H:%M:%S").timestamp()
                        info.sources["created"] = "file (TIFF DateTime tag)"
                    except ValueError:
                        pass
            info.name = output_name(fname, info.series_name, info.n_series, si)
            out.append(info)
    return out


def _page_selection(axes: str, shape: tuple, t_ax: str, channel: int, z_plane: str, max_frames):
    """Page indices [T, Z'] to read, or None when pages do not map to the axes."""
    page_axes = [a for a in axes if a not in "YXS"]
    page_shape = [shape[axes.index(a)] for a in page_axes]
    idx = np.arange(int(np.prod(page_shape)) if page_shape else 1).reshape(page_shape or [1])
    names = "".join(page_axes) or "T"
    for ax in list(names):
        if ax == t_ax:
            continue
        i = names.index(ax)
        if ax == "C":
            idx = np.take(idx, min(channel, idx.shape[i] - 1), axis=i)
        elif ax == "Z" and str(z_plane).lower() == "max":
            idx = np.moveaxis(idx, i, -1)
            names = names.replace("Z", "") + "Z"
            continue
        elif ax == "Z":
            idx = np.take(idx, min(int(z_plane or 0), idx.shape[i] - 1), axis=i)
        else:
            idx = np.take(idx, 0, axis=i)
        names = names.replace(ax, "", 1)
    if t_ax in names and names.index(t_ax) != 0:
        idx = np.moveaxis(idx, names.index(t_ax), 0)
    if t_ax not in names:
        idx = idx[None]
    if idx.ndim == 1:
        idx = idx[:, None]
    return idx[:max_frames] if max_frames else idx


def read(info: MovieInfo, opts: ImportOptions, max_frames: int | None = None) -> np.ndarray:
    import tifffile

    axes, t_ax = info.extra["axes"], info.extra["t_axis"]
    ch = info.channel_index(opts)
    with tifffile.TiffFile(info.path) as tf:
        s = tf.series[info.series]
        pages = s.pages
        shape = tuple(s.shape)
        n_pages = int(np.prod([shape[axes.index(a)] for a in axes if a not in "YXS"] or [1]))
        if len(pages) == n_pages and all(p is not None for p in pages[:1]):
            try:
                sel = _page_selection(axes, shape, t_ax, ch, opts.z_plane, max_frames)
                frames = []
                for row in sel:
                    planes = [pages[int(k)].asarray() for k in row]
                    if "S" in axes:
                        planes = [pl[..., min(ch, pl.shape[-1] - 1)] if pl.ndim == 3 else pl for pl in planes]
                    frames.append(np.max(planes, axis=0) if len(planes) > 1 else planes[0])
                return np.stack(frames)
            except Exception:
                pass   # fall back to reading the whole series
        arr = s.asarray()
    movie = select_planes(arr, axes, t_ax, ch, opts.z_plane)
    return movie[:max_frames] if max_frames else movie
