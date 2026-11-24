"""Reads sheet 'FigureS10F' (=published Supplementary Fig. 6F, normalized radial actin density)."""
from __future__ import annotations

import functools
import os
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

_M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
_XLSX = "41467_2024_46726_MOESM9_ESM.xlsx"


def _col_idx(col: str) -> int:
    n = 0
    for ch in col:
        n = n * 26 + (ord(ch) - 64)
    return n


def _f(x: object) -> float:
    return float(str(x).strip())


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


@functools.lru_cache(maxsize=2)
def _load_s10f(src_str: str):
    rows = _read_sheet(Path(src_str) / _XLSX, "FigureS10F")
    cols = sorted([c for c in rows[2] if c != "A"], key=_col_idx)
    rs = [_f(rows[2][c]) for c in cols]
    row_map = {0: 3, 10: 5, 25: 7, 50: 9, 75: 11, 100: 13, 125: 15, 150: 17}
    profiles = {conc: [_f(rows[rn][c]) for c in cols] for conc, rn in row_map.items()}
    return rs, profiles


def _interp(xs: list[float], ys: list[float], xq: float) -> float:
    if xq <= xs[0]:
        return ys[0]
    if xq >= xs[-1]:
        return ys[-1]
    for i in range(len(xs) - 1):
        if xs[i] <= xq <= xs[i + 1]:
            w = (xq - xs[i]) / (xs[i + 1] - xs[i])
            return ys[i] + w * (ys[i + 1] - ys[i])
    return ys[-1]


def run(inputs, settings):
    rs, profiles = _load_s10f(os.environ["BOXGRAPH_SOURCE"])
    conc = int(settings.get("nuc_conc_nM", 0))
    ys = profiles[conc]
    pk_i = max(range(len(ys)), key=lambda i: ys[i])
    return {
        "normalized_density": _interp(rs, ys, float(inputs["r"])),
        "peak_density": ys[pk_i],
        "peak_radius": rs[pk_i],
    }
