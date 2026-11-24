"""Reads sheet 'FigureS2B' (area-equivalent equilateral side lengths for icosphere subdivisions n=0..5)."""
from __future__ import annotations

import functools
import os
import statistics
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
def _load_s2b(src_str: str):
    rows = _read_sheet(Path(src_str) / _XLSX, "FigureS2B")
    out = {}
    for rn in range(2, 8):
        n = int(round(_f(rows[rn]["A"])))
        cols = sorted([c for c in rows[rn] if c != "A"], key=_col_idx)
        vals = [_f(rows[rn][c]) for c in cols]
        expected = 20 * (4 ** n)
        out[n] = (statistics.median(vals), min(vals), max(vals), len(vals), expected)
    return out


def run(inputs, settings):
    data = _load_s2b(os.environ["BOXGRAPH_SOURCE"])
    n = int(round(float(inputs["n"])))
    med, mn, mx, n_stored, n_exp = data[n]
    return {
        "side_median": med,
        "side_min": mn,
        "side_max": mx,
        "n_faces_stored": float(n_stored),
        "n_faces_expected": float(n_exp),
        "n_faces_missing_truncated": float(max(0, n_exp - n_stored)),
    }
