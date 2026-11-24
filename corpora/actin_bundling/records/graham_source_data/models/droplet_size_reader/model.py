"""Reads droplet diameter distributions from MOESM20 ('Fig S4e') and MOESM22 ('Fig S9c')."""
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
    ds = settings.get("dataset", "edfig4e")
    d_query = float(inputs["D"])

    if ds == "edfig4e":
        rows = _read_sheet(src / "41567_2022_1924_MOESM20_ESM.xlsx", "Fig S4e")
        vals = np.array([float(rows[r]["E"]) for r in sorted(rows) if r >= 2 and "E" in rows[r]], dtype=float)
        x0 = float(settings.get("bin_start", 0.629))
        dx = float(settings.get("bin_width", 0.3156))
        edges = [x0 + k * dx for k in range(22)]
        cnts = [int(((vals > edges[k]) & (vals <= edges[k + 1])).sum()) for k in range(21)]
        # Map digitized bar x-coordinate (bar 0 at ~1.07 um, step 0.315 um) to bin index 0..20
        k_idx = max(0, min(20, int(round((d_query - 1.07) / 0.315))))
        b_cnt = float(cnts[k_idx])
        b_freq = 100.0 * b_cnt / float(len(vals))
    elif ds.startswith("edfig9c_"):
        rows = _read_sheet(src / "41567_2022_1924_MOESM22_ESM.xlsx", "Fig S9c")
        col = {"edfig9c_3peg": "C", "edfig9c_5peg": "D", "edfig9c_10peg": "E", "edfig9c_15peg": "F"}[ds]
        vals = np.array([float(rows[r][col]) for r in sorted(rows) if r >= 3 and col in rows[r]], dtype=float)
        edges = np.arange(0.0, 4.4, 0.2)
        cnts, _ = np.histogram(vals, bins=edges)
        pcts = 100.0 * cnts / float(len(vals))
        # Published ED Fig. 9c curves plot [x, x+0.2) percentage frequency at lower bin edge x
        x_lower = edges[:-1]
        b_freq = float(np.interp(d_query, x_lower, pcts))
        b_cnt = b_freq / 100.0 * float(len(vals))
    else:
        raise ValueError(f"unknown dataset: {ds!r}")

    return {
        "n_droplets": float(len(vals)),
        "diam_median": float(np.median(vals)),
        "diam_mean": float(vals.mean()),
        "diam_sd": float(vals.std(ddof=1)),
        "diam_min": float(vals.min()),
        "diam_max": float(vals.max()),
        "bin_count": b_cnt,
        "bin_freq": b_freq,
    }
