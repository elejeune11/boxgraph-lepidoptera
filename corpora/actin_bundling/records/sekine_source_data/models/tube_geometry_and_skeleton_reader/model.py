"""Reads Source-data_Fig1j.xlsx, Source-data_Fig1k.xlsx, and Source-data_SupplementaryFig1g.xlsx from $BOXGRAPH_SOURCE."""
from __future__ import annotations

import os
import statistics
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

_M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
_TUBE_COLS = ("C", "D", "E", "F", "G", "H")
_SKEL_COLS = ("C", "D", "E", "F", "G", "H", "I", "J", "K")


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


def _find_row(rows: dict[int, dict[str, object]], t_min: float, r_lo: int, r_hi: int) -> dict[str, object] | None:
    for r in range(r_lo, r_hi + 1):
        if r in rows and "B" in rows[r] and isinstance(rows[r]["B"], (int, float)) and abs(float(rows[r]["B"]) - t_min) < 1e-6:
            return rows[r]
    return None


def run(inputs, settings):
    src = Path(os.environ["BOXGRAPH_SOURCE"])
    t_min = float(settings.get("time_min", 0))
    r1j = _read_sheet(src / "Source-data_Fig1j.xlsx", "Figure 1j")
    r1k = _read_sheet(src / "Source-data_Fig1k.xlsx", "Figure 1k")
    r1g = _read_sheet(src / "Source-data_SupplementaryFig1g.xlsx", "Supplementary Figure 1g")

    # Max absolute difference across all 36 tube cells between Fig1k and Fig1j
    diffs = [
        abs(float(r1k[r][c]) - float(r1j[r][c]))
        for r in range(7, 13)
        for c in _TUBE_COLS
    ]
    max_dup_diff = max(diffs)

    # Tube diameter / Fig1k at t_min (if present; fallback to 0 min when inspecting 10 min for skeleton)
    row_j = _find_row(r1j, t_min, 7, 12) or _find_row(r1j, 0.0, 7, 12)
    row_k = _find_row(r1k, t_min, 7, 12) or _find_row(r1k, 0.0, 7, 12)
    dj = [float(row_j[c]) for c in _TUBE_COLS]
    dk = [float(row_k[c]) for c in _TUBE_COLS]

    row_g = _find_row(r1g, t_min, 7, 17)
    if row_g is None:
        raise KeyError(f"time_min={t_min} not found in Supp. Fig. 1g")
    skels = [float(row_g[c]) for c in _SKEL_COLS if c in row_g and isinstance(row_g[c], (int, float))]

    return {
        "tube_diameter_mean": statistics.mean(dj),
        "tube_diameter_sd": statistics.stdev(dj),
        "fig1k_stored_length_mean": statistics.mean(dk),
        "fig1k_minus_fig1j_max_abs": max_dup_diff,
        "skeleton_length_mean": statistics.mean(skels),
        "skeleton_length_sd_pop": statistics.pstdev(skels),
        "skeleton_n": float(len(skels)),
    }
