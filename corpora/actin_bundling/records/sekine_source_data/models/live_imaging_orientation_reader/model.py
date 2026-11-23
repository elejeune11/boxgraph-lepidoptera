"""Reads Source-data_Fig1h.xlsx and Source-data_Fig1i.xlsx directly from $BOXGRAPH_SOURCE."""
from __future__ import annotations

import os
import statistics
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

_M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
_IMG_COLS = ("C", "D", "E", "F", "G", "H", "I", "J", "K")


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


def _row_for_time(rows: dict[int, dict[str, object]], t_min: float) -> list[float]:
    for r in range(7, 18):
        if r in rows and "B" in rows[r] and abs(float(rows[r]["B"]) - t_min) < 1e-6:
            return [float(rows[r][c]) for c in _IMG_COLS if c in rows[r] and isinstance(rows[r][c], (int, float))]
    raise KeyError(f"time_min={t_min} not found")


def run(inputs, settings):
    src = Path(os.environ["BOXGRAPH_SOURCE"])
    t_min = float(settings.get("time_min", 0))
    r1h = _read_sheet(src / "Source-data_Fig1h.xlsx", "Figure 1h")
    r1i = _read_sheet(src / "Source-data_Fig1i.xlsx", "Figure 1i")
    dirs = _row_for_time(r1h, t_min)
    cohs = _row_for_time(r1i, t_min)
    return {
        "n_samples": float(len(dirs)),
        "dominant_direction_mean": statistics.mean(dirs),
        "dominant_direction_sd_pop": statistics.pstdev(dirs),
        "coherency_mean": statistics.mean(cohs),
        "coherency_sd_pop": statistics.pstdev(cohs),
    }
