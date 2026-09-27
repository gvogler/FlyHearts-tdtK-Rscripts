"""Step 3b: combine all analyzed kymographs into the summary tables.

Port of the rest of "MAYO Screen Script No.4" and "No.5". It runs in
'<target>/balled/excellent traces' like the R script and writes the same CSV
tables (same names, columns and column order). R's .Rdata files are replaced
by CSV files (see README).
"""

from __future__ import annotations

import datetime as dt
import math
import os
import re
import shutil
from typing import Callable

import numpy as np
import pandas as pd

from . import rio
from .rstats import NA, r_mad_narm, r_mean, r_mean_narm, r_median_narm, r_sd, r_sd_narm
from .transients import METRIC_NAMES, TraceResult, analyze_trace


# ----------------------------------------------------------------------- name parsing

def _upos(s: str, n: int):
    """1-based position of the n-th underscore (None if missing)."""
    pos = [i + 1 for i, c in enumerate(s) if c == "_"]
    return pos[n - 1] if len(pos) >= n else None


def genotype_of(s: str) -> str:
    p = _upos(s, 1)
    return s[: p - 1] if p else ""


def age_of(s: str) -> str:
    p = _upos(s, 2)
    return s[p] if p and p < len(s) else ""


def sex_of(s: str) -> str:
    p = _upos(s, 3)
    return s[p - 2] if p and p >= 2 else ""


def id_of(s: str) -> str:
    p1, p2 = _upos(s, 1), _upos(s, 2)
    return s[p1:p2 - 1] if p1 and p2 else ""


def fly_part(s: str) -> str:
    """substring(x, first '_' + 1, position of 'w[fm]_' + 1)."""
    p1 = _upos(s, 1)
    m = re.search(r"w[fm]_", s)
    if not p1 or not m:
        return ""
    return s[p1: m.start() + 2]


def comment_of(s: str) -> str:
    p5 = _upos(s, 5)
    dot = s.find(".")
    if not p5 or dot < 0:
        return "none"
    return s[p5: dot]


def xpos_of(s: str) -> float:
    m = re.search(r"Xpos_(.*?).tiff", s)
    try:
        return float(m.group(1)) if m else NA
    except ValueError:
        return NA


def no_peak_number(s: str) -> str:
    """'<cxd>_peak_3_at Xpos_120...' -> '<cxd>_peak_X_at Xpos_120...'."""
    a, b = s.find("peak_"), s.find("_at")
    if a < 0 or b < 0:
        return s
    return s[: a + 5] + "X" + s[b:]


def r_sort_key(s: str):
    """Approximates R's locale sort of factor levels (case-insensitive first)."""
    return (s.lower(), s)


# ----------------------------------------------------------------------- mappings

def read_mappings(path: str) -> pd.DataFrame:
    if path.lower().endswith((".xlsx", ".xls")):
        df = pd.read_excel(path, sheet_name=0, dtype=str)
    else:
        df = pd.read_csv(path, dtype=str, encoding="utf-8-sig")
    df = df.fillna("")
    if df.shape[1] < 8:
        raise ValueError("The mappings.csv/xlsx file has missing columns")
    df.iloc[:, 0] = df.iloc[:, 0].str.replace("_", "", regex=False)
    return df


def _stock_column(crosses: pd.DataFrame):
    for c in crosses.columns:
        if re.sub(r"[^a-z]", "", str(c).lower()) == "stockcollection":
            return crosses[c]
    return None


# ----------------------------------------------------------------------- metadata table

def build_meta_data(balled: str, excellent: str) -> pd.DataFrame:
    files = rio.list_files(balled, r"_new_meta_data.csv$", recursive=True)
    rows = []
    for f in files:
        meta = pd.read_csv(os.path.join(balled, f), dtype=str, keep_default_na=False)
        rows.append(list(meta["value"].iloc[:21]) + [f])
    if not rows:
        return pd.DataFrame()
    names = list(pd.read_csv(os.path.join(balled, files[0]), dtype=str)["L1"].iloc[:21])
    md = pd.DataFrame(rows, columns=names + ["filename"])
    md["CODE"] = md["filename"].map(genotype_of)
    md["flyID"] = md["filename"].map(lambda s: s[: _upos(s, 3) - 1] if _upos(s, 3) else "")
    created = pd.to_numeric(md["created_unix_from_file"], errors="coerce")

    def day(v):
        if pd.isna(v):
            return None
        return dt.datetime.fromtimestamp(float(v), tz=dt.timezone.utc).date()

    md["recordingday"] = [day(v) for v in created]
    md["week_year"] = [d.strftime("%Y-%U") if d else None for d in md["recordingday"]]
    md["cxd_date"] = None
    md["cxd_file"] = md["filename"]
    unix_file = os.path.join(excellent, "test_unix_data.csv")
    if os.path.exists(unix_file):
        raw = pd.read_csv(unix_file, sep="\t", header=None, dtype=str).fillna("")
        raw["file"] = raw[0].map(lambda s: s.split("/")[-1])
        raw = raw.drop_duplicates(subset=0)
        raw["flyID"] = raw["file"].map(lambda s: s[: _upos(s, 3) - 1] if _upos(s, 3) else "")
        raw["V1"] = raw[0].map(lambda s: s[: s.find(" ")] if " " in s else "")
        first = raw.drop_duplicates("flyID").set_index("flyID")
        md["cxd_date"] = md["flyID"].map(first["V1"])
        md["cxd_file"] = md["flyID"].map(first["file"])
    missing = md["cxd_date"].isna()
    md.loc[missing, "cxd_date"] = md.loc[missing, "created_unix_from_file"]
    cd = pd.to_numeric(md["cxd_date"], errors="coerce")
    md["cxd_recordingday"] = [day(math.trunc(v)) if not pd.isna(v) else None for v in cd]
    md["cxd_week_year"] = [d.strftime("%U %Y") if d else None for d in md["cxd_recordingday"]]
    md["cxd_file"] = md["cxd_file"].astype(str).str.replace("_new_meta_data.csv", "cxd", regex=False)
    return md


# ----------------------------------------------------------------------- main entry

def run_analysis(excellent: str, mappings_file: str,
                 map_fn: Callable = map, log: Callable = print) -> dict:
    """Run step 3 in the 'excellent traces' folder. ``map_fn(func, iterable)``
    lets the caller run the per-kymograph analysis in parallel."""
    shutil.copy(mappings_file, os.path.join(excellent, os.path.basename(mappings_file)))
    crosses = read_mappings(mappings_file)
    c = crosses.columns
    z = pd.DataFrame({"CODE": crosses[c[0]], "cross": crosses[c[1]], "type": crosses[c[2]],
                      "fly": crosses[c[5]], "human": crosses[c[6]]})
    zfirst = z.drop_duplicates("CODE").set_index("CODE")

    jpgs = rio.list_files(excellent, r"\.jpg$")
    if not jpgs:
        raise RuntimeError("No traced kymographs (*.jpg) in the 'excellent traces' folder")
    fl = pd.DataFrame({"jpg": jpgs})
    fl["file"] = fl["jpg"].str.replace(r"_traced.jpg$", ".csv", regex=True)
    fl["Genotype_filename"] = fl["file"].map(genotype_of)
    fl["Age_filename"] = fl["file"].map(age_of)
    fl["Sex_filename"] = fl["file"].map(sex_of)
    fl["comment"] = fl["file"].map(comment_of)
    fl["index"] = np.arange(1, len(fl) + 1)
    fl["CODED"] = fl["Genotype_filename"] + fl["Age_filename"] + fl["Sex_filename"]
    fl["CODE"] = fl["Genotype_filename"]
    for col, src in (("cross", "cross"), ("type", "type"), ("symbol_fly", "fly"), ("symbol_human", "human")):
        fl[col] = fl["CODE"].map(zfirst[src])
    fl["Genotype"] = fl["CODED"]
    fl["IDnumber"] = pd.to_numeric(fl["file"].map(id_of), errors="coerce")
    fl["flyID"] = fl["Genotype_filename"] + "_" + fl["file"].map(fly_part)

    # direction files live one folder up ('balled')
    balled = os.path.dirname(os.path.normpath(excellent))
    dfiles = rio.list_files(balled, r"_directionmarks.csv")
    dmap = {rio.cxd_stem(f): os.path.join(balled, f) for f in dfiles}

    # ---- analyze every kymograph (parallel) ------------------------------
    tasks = [(excellent, f, dmap.get(rio.cxd_stem(f))) for f in fl["file"]]
    results: list[TraceResult] = list(map_fn(_analyze_task, tasks))
    for r in results:
        if not r.ok:
            log(f"  skipped {r.csv}: {r.message}")
    ok = np.array([r.ok for r in results])
    res_by_index = {int(ix): r for ix, r in zip(fl["index"], results)}

    # ---- genotype groups (only analyzable kymographs) --------------------
    fl_ok = fl[ok]
    if fl_ok.empty:
        raise RuntimeError("None of the kymographs in 'excellent traces' could be analyzed (see the log)")
    codes = sorted(fl_ok["Genotype"].unique(), key=r_sort_key)
    groups = [(i + 1, code, list(fl_ok.loc[fl_ok["Genotype"] == code, "index"])) for i, code in enumerate(codes)]

    intervals, transients, final_tr = [], [], []
    for i, code, idxs in groups:
        for j, ix in enumerate(idxs, start=1):
            r = res_by_index[ix]
            iv = r.intervals.copy()
            iv["i"], iv["j"], iv["Index"] = i, j, ix
            intervals.append(iv)
            tr = r.transients.copy()
            tr["i"], tr["j"], tr["index"] = i, j, ix
            transients.append(tr[["time", "distance", "diameter", "timestamp", "delta", "from_EDD",
                                  "velocity", "i", "j", "index", "beat"]])
            for m in r.metrics:
                final_tr.append([i, j, ix] + list(m))
    intervals_final = pd.concat(intervals, ignore_index=True) if intervals else pd.DataFrame()
    transients_final = pd.concat(transients, ignore_index=True) if transients else pd.DataFrame()
    final_transients = pd.DataFrame(final_tr, columns=["i", "j", "index"] + METRIC_NAMES)

    rio.write_r_csv(intervals_final, os.path.join(excellent, "intervals_final.csv"))
    transients_final.to_csv(os.path.join(excellent, "transients_final.csv.gz"), index=False, na_rep="NA")
    rio.write_r_csv(final_transients, os.path.join(excellent, "final_transients.csv"))

    # ---- per-fly values ----------------------------------------------------
    per_file = pd.DataFrame({
        "index": fl["index"],
        "total_time": [res_by_index[ix].total_time for ix in fl["index"]],
        "time_used": [res_by_index[ix].time_used for ix in fl["index"]],
        "reversals": [res_by_index[ix].reversals for ix in fl["index"]],
        "anterior_beat_percent": [res_by_index[ix].anterior_beat_percent for ix in fl["index"]],
        "speed_anterograde": [res_by_index[ix].speed_anterograde for ix in fl["index"]],
        "speed_retrograde": [res_by_index[ix].speed_retrograde for ix in fl["index"]],
    }).set_index("index")

    def fly_values(ix):
        r = res_by_index[ix]
        tt = round(r.total_time, 2) if not math.isnan(r.total_time) else NA
        return {
            "HR": r.n_identicals / tt if tt else NA,
            "EDD": -r.edd, "ESD": -r.esd, "FS": r.fs,
            "min_velocity": r_mean_narm(r.max_neg_velocity),
            "max_velocity": r_mean_narm(r.max_velocity),
            "coverage": r.time_used / r.total_time if r.total_time else NA,
        }

    master_rows, long_rows = [], []
    for i, code, idxs in groups:
        vals = [fly_values(ix) for ix in idxs]
        orig = fl.loc[fl["Genotype"] == code, "CODE"].iloc[0]
        row = {"Code": code, "n": len(idxs),
               "Cross": zfirst["cross"].get(orig, NA), "type": zfirst["type"].get(orig, NA),
               "human": zfirst["human"].get(orig, NA), "fly": zfirst["fly"].get(orig, NA),
               "AgeSex": code[-2:], "Sex": code[-1:], "Age": code[-2:][:1]}
        cov = [v["coverage"] for v in vals]
        row["coverage"], row["SDcoverage"] = r_mean(cov), r_sd(cov)
        for key, mean_name, sd_name in (("HR", "HR", "SDHR"), ("EDD", "EDD", "SDEDD"), ("ESD", "ESD", "SDESD"),
                                        ("FS", "FS", "SDFS"), ("max_velocity", "max_velocity", "SDMV"),
                                        ("min_velocity", "min_velocity", "SDMinV")):
            xs = [v[key] for v in vals]
            row[mean_name], row[sd_name] = r_mean(xs), r_sd(xs)
        master_rows.append(row)
        for ix, v in zip(idxs, vals):
            long_rows.append({"index": ix, "L1": code, **v})
    master = pd.DataFrame(master_rows)
    rio.write_r_csv(master, os.path.join(excellent, "master_table.csv"))

    long = pd.DataFrame(long_rows)
    long["CODE"] = long["L1"].map(lambda s: s[:-2] if s[-1:] in ("m", "f") else "")
    long["human"] = long["CODE"].map(zfirst["human"])
    long["fly"] = long["CODE"].map(zfirst["fly"])
    long["type"] = long["CODE"].map(zfirst["type"])
    for col, fname, label in (("EDD", "EDD_table.csv", "EDD"), ("ESD", "ESD_table.csv", "ESD"),
                              ("FS", "FS_table.csv", "FS"), ("HR", "HR_table.csv", "HR"),
                              ("min_velocity", "min_velocity_data.csv", "min_velocity"),
                              ("max_velocity", "max_velocity_table.csv", "max.velocity")):
        t = long[[col, "L1", "CODE", "human", "fly", "type"]].rename(columns={col: label})
        rio.write_r_csv(t, os.path.join(excellent, fname))

    # ---- metadata ------------------------------------------------------------
    meta = build_meta_data(balled, excellent)
    if not meta.empty:
        rio.write_r_csv(meta, os.path.join(excellent, "meta_data_all.csv"))

    # ---- all_data_table ------------------------------------------------------
    jpg_of = fl.set_index("index")["jpg"]
    ad = pd.DataFrame({
        "file": [jpg_of[ix] for ix in long["index"]],
        "EDD": long["EDD"], "ESD": long["ESD"], "FS": long["FS"], "HR": long["HR"],
        "min.velocity": long["min_velocity"], "max.velocity": long["max_velocity"],
        "CODE_long": long["L1"], "CODE": long["CODE"], "human_gene": long["human"],
        "fly_gene": long["fly"], "type": long["type"],
    })
    ad["cross"] = ad["CODE"].map(zfirst["cross"])
    ad["flyID"] = ad["CODE"] + "_" + ad["file"].map(fly_part)
    ad["ID_number"] = pd.to_numeric(ad["file"].map(id_of), errors="coerce")
    ad["Xpos"] = ad["file"].map(xpos_of)
    ad["cxd_file"] = ad["file"].map(rio.cxd_stem)
    ad["_order"] = long["index"].to_numpy()
    if not meta.empty:
        m2 = meta.drop(columns=["CODE", "flyID"])
        for col in m2.columns[:21]:
            conv = pd.to_numeric(m2[col], errors="coerce")
            if conv.notna().all():
                m2[col] = conv
        ad = ad.merge(m2, on="cxd_file", how="left")
    else:
        for col in ["sizeX", "sizeY", "sizeZ", "sizeC", "sizeT", "pixelType", "bitsPerPixel", "imageCount",
                    "dimensionOrder", "orderCertain", "rgb", "littleEndian", "interleaved", "falseColor",
                    "metadataComplete", "thumbnail", "series", "resolutionLevel", "time_interval",
                    "resolution", "created_unix_from_file", "filename", "recordingday", "week_year",
                    "cxd_date", "cxd_recordingday", "cxd_week_year"]:
            ad[col] = NA
    ad = ad.sort_values("cxd_file", kind="stable").reset_index(drop=True)   # merge() sorts by key
    ad = ad.drop(columns=["_order"])
    ad["CODE.x"] = ad["CODE"]
    ad["CODE.y"] = ad["CODE"]
    meta_cols = ["sizeX", "sizeY", "sizeZ", "sizeC", "sizeT", "pixelType", "bitsPerPixel", "imageCount",
                 "dimensionOrder", "orderCertain", "rgb", "littleEndian", "interleaved", "falseColor",
                 "metadataComplete", "thumbnail", "series", "resolutionLevel", "time_interval",
                 "resolution", "created_unix_from_file", "filename"]
    all_data = ad[["flyID", "file", "EDD", "ESD", "FS", "HR", "min.velocity", "max.velocity", "CODE_long",
                   "CODE.x", "human_gene", "fly_gene", "type", "cross", "ID_number", "Xpos"] + meta_cols
                  + ["CODE.y", "recordingday", "week_year", "cxd_date", "cxd_file", "cxd_recordingday",
                     "cxd_week_year"]].copy()

    # AI and MAD of tt10r per kymograph
    ai = (intervals_final.groupby(["i", "j", "Index"])["HP"]
          .agg(sd=r_sd_narm, median=r_median_narm).reset_index())
    ai["AI"] = ai["sd"] / ai["median"]
    mad = final_transients.groupby(["i", "j", "index"])["tt10r"].agg(r_mad_narm).reset_index()
    mad = mad.rename(columns={"tt10r": "MAD_tt10r"})
    ai = ai.merge(mad, left_on=["i", "j", "Index"], right_on=["i", "j", "index"], how="left")
    ai_by_index = ai.set_index("Index")

    sv = all_data[["flyID", "file", "EDD", "ESD", "HR", "Xpos"]].copy()
    sv["SV"] = math.pi / 4 * (sv["EDD"] ** 2 - sv["ESD"] ** 2)
    sv["CO"] = sv["SV"] * sv["HR"]
    thr_x = math.floor(np.mean(pd.unique(sv["Xpos"].dropna()))) if sv["Xpos"].notna().any() else NA
    sv["segment"] = np.where(sv["Xpos"] > thr_x, "posterior", "anterior")
    sv_by_file = sv.drop_duplicates("file").set_index("file")

    npar = pd.DataFrame({"jpg": fl["jpg"], "Genotype": fl["Genotype"], "Genotype_filename": fl["Genotype_filename"]})
    npar["AI"] = fl["index"].map(ai_by_index["AI"]).to_numpy()
    npar["MAD_tt10r"] = fl["index"].map(ai_by_index["MAD_tt10r"]).to_numpy()
    for col in ("SV", "CO", "segment"):
        npar[col] = npar["jpg"].map(sv_by_file[col])
    npar["jpg_noX"] = npar["jpg"].map(no_peak_number)
    npar["anterograde_percent"] = fl["index"].map(per_file["anterior_beat_percent"]).to_numpy()
    npar["anterograde_speed"] = fl["index"].map(per_file["speed_anterograde"]).to_numpy()
    npar["retrograde_speed"] = fl["index"].map(per_file["speed_retrograde"]).to_numpy()
    all_data_added = all_data.merge(npar, left_on="file", right_on="jpg", how="left").drop(columns=["jpg"])
    rio.write_r_csv(all_data_added, os.path.join(excellent, "all_data_table.csv"))

    # ---- per-kymograph transient summary --------------------------------------
    g = final_transients.groupby(["i", "j", "index"])[METRIC_NAMES]
    gpf = pd.concat([g.agg(r_median_narm).add_suffix("_median"), g.agg(r_sd_narm).add_suffix("_sd")], axis=1)
    gpf = gpf.reset_index()
    gpf["CODE"] = gpf["index"].map(fl.set_index("index")["Genotype_filename"])
    gpf["jpg"] = gpf["index"].map(jpg_of)
    icols = ["starts", "peaks", "ends", "SI", "deltaSI", "deltaSIplus", "HP", "DI"]
    gi = intervals_final.groupby(["i", "j", "Index"])[icols]
    ipf = pd.concat([gi.agg(r_median_narm).add_suffix("_median"), gi.agg(r_sd_narm).add_suffix("_sd")], axis=1)
    ipf = ipf.reset_index()
    mads = intervals_final.groupby(["i", "j", "Index"])[["SI", "HP", "DI"]].agg(r_mad_narm)
    mads.columns = ["MAD_SI", "MAD_HP", "MAD_DI"]
    mads = mads.reset_index()
    mt = gpf.merge(ipf, left_on=["i", "j", "index"], right_on=["i", "j", "Index"], how="left").drop(columns=["Index"])
    mt["relaxtime"] = mt["HP_median"] - mt["peaked_median"]
    mt = mt.merge(mads, left_on=["i", "j", "index"], right_on=["i", "j", "Index"], how="left").drop(columns=["Index"])

    fad = all_data_added.merge(mt.drop(columns=["i", "j", "index"]), left_on="file", right_on="jpg", how="left")
    fad = fad.drop(columns=["jpg"]).drop_duplicates()
    fad["Age"] = fad["file"].map(age_of)
    fad["Sex"] = fad["file"].map(sex_of)
    fad = fad.drop_duplicates("file", keep="first")
    fad["file"] = fad["file"].map(no_peak_number)
    order = ["CODE", "CODE_long", "flyID", "type", "cross", "cxd_file", "cxd_recordingday", "cxd_week_year",
             "human_gene", "fly_gene", "Xpos", "Age", "Sex", "EDD", "ESD", "FS", "HR", "HP_median",
             "DI_median", "relaxtime", "SI_median", "peaked_median", "min.velocity", "max.velocity",
             "MAD_SI", "MAD_HP", "MAD_DI"] + [m + "_median" for m in METRIC_NAMES[1:]] + [
             "deltaSI_median", "deltaSIplus_median", "AI", "MAD_tt10r", "SV", "CO", "anterograde_percent",
             "segment", "file", "ID_number", "time_interval", "resolution", "filename", "cxd_date",
             "Genotype", "Genotype_filename", "jpg_noX", "anterograde_speed", "retrograde_speed"]
    fad = fad[order].copy()
    stock = _stock_column(crosses)
    if stock is not None:
        smap = pd.Series(stock.values, index=crosses[c[0]]).groupby(level=0).first()
        fad["stocks"] = fad["CODE"].map(smap)
    else:
        fad["stocks"] = "empty"
    fad = fad.drop_duplicates().reset_index(drop=True)
    rio.write_r_csv(fad, os.path.join(excellent, "final_all_data.csv"))

    # ---- per fly ----------------------------------------------------------------
    med_cols = order[13:44] + ["anterograde_speed", "retrograde_speed"]
    id_cols = order[0:10] + ["Age", "Sex", "time_interval", "resolution", "filename", "cxd_date",
                             "Genotype", "Genotype_filename"]

    def per_fly(df: pd.DataFrame, extra: str) -> pd.DataFrame:
        med = df.groupby("flyID", sort=True)[med_cols].agg(r_median_narm).reset_index()
        ids = df[[cc for cc in id_cols] + [extra]].drop(columns=["flyID"])
        ids.insert(0, "flyID", df["flyID"])
        out = med.merge(ids, on="flyID", how="left", suffixes=(".x", ".y"))
        return out.drop_duplicates().reset_index(drop=True)

    pf = per_fly(fad, "stocks")
    rev = fl[["flyID"]].copy()
    rev["reversals"] = fl["index"].map(per_file["reversals"]).to_numpy()
    rev = rev[ok].drop_duplicates()
    rev = rev.groupby("flyID")["reversals"].agg(reversals_median=r_median_narm, reversals_SD=r_sd_narm).reset_index()
    pf = pf.merge(rev, on="flyID", how="left")
    rio.write_r_csv(pf, os.path.join(excellent, "final_all_data_per_fly.csv"))

    for seg in ("anterior", "posterior"):
        part = fad[fad["segment"] == seg]
        rio.write_r_csv(part, os.path.join(excellent, f"final_all_data_{seg}.csv"))
        # (the R script joins column 55 = anterograde_speed here instead of 'stocks')
        rio.write_r_csv(per_fly(part, "anterograde_speed"),
                        os.path.join(excellent, f"final_all_data_per_fly_{seg}.csv"))

    n_ok = int(ok.sum())
    return {"kymographs": len(fl), "analyzed": n_ok, "genotypes": len(groups),
            "skipped": [(r.csv, r.message) for r in results if not r.ok]}


def _analyze_task(args) -> TraceResult:
    folder, csv_name, dfile = args
    return analyze_trace(folder, csv_name, dfile)
