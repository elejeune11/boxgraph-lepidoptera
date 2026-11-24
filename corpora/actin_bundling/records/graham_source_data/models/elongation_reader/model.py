"""Reads linear droplet-bundle elongation trajectories from MOESM16 ('Fig 5h')."""
from __future__ import annotations

import os
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np

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
        rows: dict[int, dict[str, object]] = {}
        for c in ET.fromstring(z.read(target)).iter(_M + "c"):
            ref = c.get("r")
            col = "".join(ch for ch in ref if ch.isalpha())
            rn = int("".join(ch for ch in ref if ch.isdigit()))
            kind, v = c.get("t"), c.find(_M + "v")
            if v is not None and kind is None:
                rows.setdefault(rn, {})[col] = float(v.text)
    return rows


def run(inputs, settings):
    src = Path(os.environ["BOXGRAPH_SOURCE"])
    rows = _read_sheet(src / "41567_2022_1924_MOESM16_ESM.xlsx", "Fig 5h")
    t_query = float(inputs["t"])

    dcols = ["B", "C", "D", "E", "F", "H", "I", "J", "K", "L", "N", "O", "P", "Q", "R"]
    trial_groups = [dcols[0:5], dcols[5:10], dcols[10:15]]
    data_rows = [r for r in sorted(rows) if r >= 4 and "A" in rows[r]]
    ts = np.array([float(rows[r]["A"]) for r in data_rows], dtype=float)
    avg_y = np.array([np.mean([float(rows[r][c]) for c in dcols if c in rows[r]]) for r in data_rows], dtype=float)

    # Through-origin regression y = m * t on average curve
    m_orig = float(np.sum(ts * avg_y) / np.sum(ts ** 2))
    r2_orig = float(1.0 - np.sum((avg_y - m_orig * ts) ** 2) / np.sum((avg_y - avg_y.mean()) ** 2))

    # Free-intercept OLS regression y = m * t + b on average curve
    m_ols, b_ols = np.polyfit(ts, avg_y, 1)
    m_ols, b_ols = float(m_ols), float(b_ols)
    r2_ols = float(1.0 - np.sum((avg_y - (m_ols * ts + b_ols)) ** 2) / np.sum((avg_y - avg_y.mean()) ** 2))

    # Individual droplet slopes
    sl_orig, sl_ols = [], []
    for c in dcols:
        pts = np.array([(float(rows[r]["A"]), float(rows[r][c])) for r in data_rows if c in rows[r]], dtype=float)
        sl_orig.append(float(np.sum(pts[:, 0] * pts[:, 1]) / np.sum(pts[:, 0] ** 2)) * 1000.0)
        sl_ols.append(float(np.polyfit(pts[:, 0], pts[:, 1], 1)[0]) * 1000.0)

    # 3-trial average through-origin slopes
    sl_trial_orig = []
    for grp in trial_groups:
        ty = np.array([np.mean([float(rows[r][c]) for c in grp if c in rows[r]]) for r in data_rows], dtype=float)
        sl_trial_orig.append(float(np.sum(ts * ty) / np.sum(ts ** 2)) * 1000.0)

    end_rates = [float(rows[data_rows[-1]][c]) / float(ts[-1]) * 1000.0 for c in dcols]
    r_20 = next(r for r in data_rows if abs(float(rows[r]["A"]) - 20.0) < 1e-6)
    n_at_20 = sum(1 for c in dcols if c in rows[r_20])

    return {
        "growth_origin_fit": m_orig * t_query,
        "growth_ols_fit": m_ols * t_query + b_ols,
        "growth_raw_avg": float(np.interp(t_query, ts, avg_y)),
        "slope_avg_origin": m_orig,
        "r2_avg_origin": r2_orig,
        "slope_avg_ols": m_ols,
        "r2_avg_ols": r2_ols,
        "slope_ind_origin_mean": float(np.mean(sl_orig)),
        "slope_ind_origin_sd": float(np.std(sl_orig, ddof=1)),
        "slope_ind_ols_mean": float(np.mean(sl_ols)),
        "slope_ind_ols_sd": float(np.std(sl_ols, ddof=1)),
        "slope_3trial_origin_sd_pop": float(np.std(sl_trial_orig, ddof=0)),
        "end_rate_mean": float(np.mean(end_rates)),
        "end_rate_sd": float(np.std(end_rates, ddof=1)),
        "n_droplets": float(len(dcols)),
        "n_at_20s": float(n_at_20),
    }
