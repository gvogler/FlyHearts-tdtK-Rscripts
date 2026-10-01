"""Manual review of traced kymographs (M-modes).

The automatic quality control puts traces into 'excellent', 'good' and 'bad
traces'; only 'excellent traces' are analyzed in step 3. Some traces that the
automatic scoring rejects are still usable (and occasionally an 'excellent'
one is not). The user's decisions are stored in 'balled/manual_curation.csv':

    csv, decision, automatic, changed
    <...>.tiff.csv, include | exclude, excellent | rescued | good | bad, <date>

'include' copies the trace (CSV + traced JPG) into 'excellent traces',
'exclude' removes it from there. The decisions are re-applied every time the
quality control runs again, so they are never lost. Step 3 must be re-run
after changes.
"""

from __future__ import annotations

import datetime as dt
import os
import shutil

import pandas as pd

from . import rio

CURATION = "manual_curation.csv"
EXCELLENT = "excellent traces"
AUTO_LABEL = {"excellent": "excellent", "good": "good", "n.t.": "bad"}


def automatic_selection(qc: pd.DataFrame) -> dict:
    """csv -> automatic class: 'excellent', 'rescued' (good trace of a movie without any
    excellent trace, copied to 'excellent traces' by the QC), 'good' or 'bad'."""
    if qc is None or qc.empty:
        return {}
    q = qc.copy()
    q["movie"] = q["jpeg"].map(rio.movie_stem)
    counts = q.groupby(["movie", "quality"]).size().unstack(fill_value=0)
    rescued_movies = set()
    if "excellent" in counts and "good" in counts:
        rescued_movies = set(counts.index[(counts["excellent"] == 0) & (counts["good"] > 0)])
    out = {}
    for _, r in q.iterrows():
        a = AUTO_LABEL.get(r["quality"], "bad")
        if a == "good" and r["movie"] in rescued_movies:
            a = "rescued"
        out[r["csv"]] = a
    return out


def load_decisions(balled: str) -> dict:
    path = os.path.join(balled, CURATION)
    if not os.path.exists(path):
        return {}
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    return {r["csv"]: r["decision"] for _, r in df.iterrows() if r["decision"] in ("include", "exclude")}


def _save_decisions(balled: str, decisions: dict, auto: dict) -> None:
    path = os.path.join(balled, CURATION)
    old = {}
    if os.path.exists(path):
        prev = pd.read_csv(path, dtype=str, keep_default_na=False)
        old = {r["csv"]: r.get("changed", "") for _, r in prev.iterrows()}
    now = dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    rows = [{"csv": c, "decision": d, "automatic": auto.get(c, ""), "changed": old.get(c) or now}
            for c, d in sorted(decisions.items())]
    pd.DataFrame(rows, columns=["csv", "decision", "automatic", "changed"]).to_csv(path, index=False)


def _jpeg(csv_name: str) -> str:
    return csv_name[: -len(".csv")] + "_traced.jpg" if csv_name.endswith(".csv") else csv_name


def _set_in_excellent(balled: str, csv_name: str, present: bool) -> None:
    target = os.path.join(balled, EXCELLENT)
    os.makedirs(target, exist_ok=True)
    for f in (csv_name, _jpeg(csv_name)):
        dst = os.path.join(target, f)
        if present:
            src = os.path.join(balled, f)
            if os.path.exists(src) and not os.path.exists(dst):
                shutil.copy(src, dst)
        elif os.path.exists(dst):
            os.remove(dst)


def _read_qc(balled: str) -> pd.DataFrame:
    path = os.path.join(balled, "Quality_control.csv")
    return rio.read_r_csv(path) if os.path.exists(path) else pd.DataFrame()


def apply_decisions(balled: str) -> int:
    """Make 'excellent traces' follow the manual decisions; returns the number of decisions."""
    decisions = load_decisions(balled)
    for c, d in decisions.items():
        _set_in_excellent(balled, c, d == "include")
    return len(decisions)


def set_decision(balled: str, csv_names, decision: str) -> None:
    """decision: 'include', 'exclude' or 'auto' (forget the manual decision)."""
    if decision not in ("include", "exclude", "auto"):
        raise ValueError(decision)
    auto = automatic_selection(_read_qc(balled))
    decisions = load_decisions(balled)
    for c in csv_names:
        if decision == "auto":
            decisions.pop(c, None)
            _set_in_excellent(balled, c, auto.get(c) in ("excellent", "rescued"))
        else:
            decisions[c] = decision
            _set_in_excellent(balled, c, decision == "include")
    _save_decisions(balled, decisions, auto)


def trace_table(balled: str) -> pd.DataFrame:
    """All traced kymographs with automatic class, manual decision and QC numbers."""
    qc = _read_qc(balled)
    cols = ["csv", "jpeg", "movie", "Xpos", "automatic", "decision", "in_analysis",
            "sd", "autocorr_up", "autocorr_down", "dist_cor"]
    if qc.empty:
        return pd.DataFrame(columns=cols)
    auto = automatic_selection(qc)
    decisions = load_decisions(balled)
    exc = os.path.join(balled, EXCELLENT)
    t = qc[["csv", "jpeg", "Xpos", "sd", "autocorr_up", "autocorr_down", "dist_cor"]].copy()
    t["movie"] = t["jpeg"].map(rio.movie_stem)
    t["automatic"] = t["csv"].map(auto)
    t["decision"] = t["csv"].map(lambda c: decisions.get(c, ""))
    t["in_analysis"] = t["csv"].map(lambda c: os.path.exists(os.path.join(exc, c)))
    return t[cols].reset_index(drop=True)


def analysis_is_stale(balled: str) -> bool:
    """True when traces were added/removed after the last beat analysis (step 3)."""
    exc = os.path.join(balled, EXCELLENT)
    result = os.path.join(exc, "final_all_data.csv")
    cur = os.path.join(balled, CURATION)
    if not os.path.exists(cur):
        return False
    return not os.path.exists(result) or os.path.getmtime(cur) > os.path.getmtime(result)
