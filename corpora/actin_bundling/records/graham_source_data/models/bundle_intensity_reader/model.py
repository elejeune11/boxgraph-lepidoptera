"""Reads normalized fluorescence intensity trajectories from MOESM17 ('Fig 6f', 'Fig 6g')."""
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
    series = settings.get("series", "fig6f_d1_actin")
    add_one = bool(settings.get("add_one", False))
    snap_nearest = bool(settings.get("snap_nearest", False))
    t_query = float(inputs["t"])

    wb_path = src / "41567_2022_1924_MOESM17_ESM.xlsx"
    if series.startswith("fig6f_"):
        rows = _read_sheet(wb_path, "Fig 6f")
        col_map = {
            "fig6f_d1_actin": "C",
            "fig6f_d1_vasp": "D",
            "fig6f_d2_actin": "F",
            "fig6f_d2_vasp": "G",
        }
        col = col_map[series]
        data_rows = [r for r in sorted(rows) if r >= 3 and "B" in rows[r] and col in rows[r]]
        ts = np.array([float(rows[r]["B"]) for r in data_rows], dtype=float)
        ys = np.array([float(rows[r][col]) for r in data_rows], dtype=float)
        if snap_nearest:
            idx = int(np.argmin(np.abs(ts - t_query)))
            val = float(ys[idx])
        else:
            val = float(np.interp(t_query, ts, ys))
        n_drop = 1.0
    elif series == "fig6g_avg":
        rows = _read_sheet(wb_path, "Fig 6g")
        gcols = ["B", "C", "D", "E", "F", "H", "I", "J", "K", "L", "N", "O", "P", "Q", "R"]
        data_rows = [r for r in sorted(rows) if r >= 4 and "A" in rows[r]]
        ts = np.array([float(rows[r]["A"]) for r in data_rows], dtype=float)
        ys = np.array([np.mean([float(rows[r][c]) for c in gcols if c in rows[r]]) for r in data_rows], dtype=float)
        val = float(np.interp(t_query, ts, ys))
        if add_one:
            val += 1.0
        n_drop = float(len(gcols))
    else:
        raise ValueError(f"unknown series: {series!r}")

    return {"normalized_intensity": val, "n_droplets": n_drop}
