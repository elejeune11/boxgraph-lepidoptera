"""Reads FRAP recovery traces from MOESM12 ('Fig 1f') and MOESM17 ('Fig 6b') directly from $BOXGRAPH_SOURCE."""
from __future__ import annotations

import math
import os
import statistics
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
from scipy.optimize import curve_fit

_M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def _read_sheet(path: Path, sheet_name: str) -> dict[int, dict[str, object]]:
    with zipfile.ZipFile(path) as z:
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        relmap = {r.get("Id"): r.get("Target") for r in rels}
        rid = next(s.get(_R + "id") for s in wb.iter(_M + "sheet") if s.get("name") == sheet_name)
        target = relmap[rid].lstrip("/")
        if not target.startswith("xl/"):
            target = "xl/" + target
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            for si in ET.fromstring(z.read("xl/sharedStrings.xml")).iter(_M + "si"):
                shared.append("".join(t.text or "" for t in si.iter(_M + "t")))
        rows: dict[int, dict[str, object]] = {}
        for c in ET.fromstring(z.read(target)).iter(_M + "c"):
            ref = c.get("r")
            col = "".join(ch for ch in ref if ch.isalpha())
            rn = int("".join(ch for ch in ref if ch.isdigit()))
            kind, v = c.get("t"), c.find(_M + "v")
            if kind == "inlineStr":
                val = "".join(t.text or "" for t in c.iter(_M + "t"))
            elif v is None:
                val = None
            elif kind == "s":
                val = shared[int(v.text)]
            elif kind in ("str", "e"):
                val = v.text
            elif kind == "b":
                val = v.text == "1"
            else:
                val = float(v.text)
            if val is not None:
                rows.setdefault(rn, {})[col] = val
    return rows


def run(inputs, settings):
    src = Path(os.environ["BOXGRAPH_SOURCE"])
    ds = settings.get("dataset", "fig1f")
    if ds == "fig1f":
        rows = _read_sheet(src / "41567_2022_1924_MOESM12_ESM.xlsx", "Fig 1f")
        tcol, cols, r0, r_bleach = "B", ["C", "D", "E", "F", "G", "H", "I", "J", "K"], 4, 7
    elif ds == "fig6b_3peg":
        rows = _read_sheet(src / "41567_2022_1924_MOESM17_ESM.xlsx", "Fig 6b")
        tcol, cols, r0, r_bleach = "A", ["B", "C", "D", "E", "F"], 5, 12
    elif ds == "fig6b_10peg":
        rows = _read_sheet(src / "41567_2022_1924_MOESM17_ESM.xlsx", "Fig 6b")
        tcol, cols, r0, r_bleach = "H", ["I", "J", "K", "L", "M"], 5, 12
    else:
        raise ValueError(f"unknown dataset: {ds!r}")

    data_rows = [r for r in sorted(rows) if r >= r0 and tcol in rows[r]]
    last_row = data_rows[-1]
    final_vals = [float(rows[last_row][c]) * 100.0 for c in cols if c in rows[last_row]]
    all_last = []
    for c in cols:
        rs = [r for r in data_rows if c in rows[r]]
        all_last.append(float(rows[rs[-1]][c]) * 100.0)

    # Interpolate average at t = 120 s
    ts_avg, ys_avg = [], []
    for r in data_rows:
        vs = [float(rows[r][c]) * 100.0 for c in cols if c in rows[r]]
        if vs:
            ts_avg.append(float(rows[r][tcol]))
            ys_avg.append(statistics.mean(vs))
    rec_120 = float(np.interp(120.0, ts_avg, ys_avg))

    # Half-times for Fig 1f (single-exponential fit A*(1-exp(-k*t)) and 50% crossing)
    thalf_exp, thalf_50 = [], []
    if ds == "fig1f":
        for c in cols:
            pts = [(float(rows[r][tcol]), float(rows[r][c])) for r in data_rows if r >= r_bleach and c in rows[r]]
            ts = np.array([p[0] for p in pts])
            ys = np.array([p[1] for p in pts])
            popt, _ = curve_fit(lambda t, a, k: a * (1.0 - np.exp(-k * t)), ts, ys, p0=[0.95, 0.04])
            thalf_exp.append(math.log(2.0) / float(popt[1]))
            for (t0, y0), (t1, y1) in zip(pts, pts[1:]):
                if y0 <= 0.5 <= y1:
                    thalf_50.append(t0 + (0.5 - y0) / (y1 - y0) * (t1 - t0))
                    break
    else:
        thalf_exp = [0.0, 0.0]
        thalf_50 = [0.0, 0.0]

    sd_fin = statistics.stdev(final_vals)
    return {
        "n_trials_total": float(len(cols)),
        "n_trials_final": float(len(final_vals)),
        "final_recovery_mean": statistics.mean(final_vals),
        "final_recovery_sd": sd_fin,
        "final_recovery_sem": sd_fin / math.sqrt(len(final_vals)),
        "all_trials_last_mean": statistics.mean(all_last),
        "all_trials_last_sd": statistics.stdev(all_last),
        "recovery_120s_mean": rec_120,
        "t_half_exp_mean": statistics.mean(thalf_exp),
        "t_half_exp_sd": statistics.stdev(thalf_exp),
        "t_half_halfmax_mean": statistics.mean(thalf_50),
        "t_half_halfmax_sd": statistics.stdev(thalf_50),
    }
