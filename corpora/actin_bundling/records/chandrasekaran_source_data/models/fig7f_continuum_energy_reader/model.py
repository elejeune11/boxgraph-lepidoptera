"""Reads sheet 'Figure7F' (2D continuum bending and surface energy vs aspect ratio)."""
from __future__ import annotations

import functools
import math
import os
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

_M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
_XLSX = "41467_2024_46726_MOESM9_ESM.xlsx"


def _f(x: object) -> float:
    return float(str(x).strip())


def _ellipe_agm(m: float) -> float:
    """Complete elliptic integral of the second kind E(m) via AGM (pure stdlib)."""
    a = 1.0
    b = math.sqrt(max(0.0, 1.0 - m))
    c_sq_sum = 0.5 * m
    p = 1.0
    for _ in range(30):
        an = 0.5 * (a + b)
        bn = math.sqrt(a * b)
        cn = 0.5 * (a - b)
        c_sq_sum += p * cn * cn
        a, b = an, bn
        p *= 2.0
        if cn < 1e-15:
            break
    return (math.pi / (2.0 * a)) * (1.0 - c_sq_sum)


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
def _load_7f(src_str: str):
    rows = _read_sheet(Path(src_str) / _XLSX, "Figure7F")
    xs, eb, es = [], [], []
    max_perim_diff = 0.0
    for rn in range(2, 21):
        ar = _f(rows[rn]["A"])
        b_v = _f(rows[rn]["B"])
        s_v = _f(rows[rn]["C"])
        xs.append(ar)
        eb.append(b_v)
        es.append(s_v)
        perim = 4.0 * math.sqrt(ar) * _ellipe_agm(1.0 - 1.0 / (ar * ar))
        max_perim_diff = max(max_perim_diff, abs(s_v - perim))
    return xs, eb, es, max_perim_diff


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
    xs, eb, es, max_perim_diff = _load_7f(os.environ["BOXGRAPH_SOURCE"])
    ar = float(inputs["aspect_ratio"])
    return {
        "bending_energy": _interp(xs, eb, ar),
        "surface_energy": _interp(xs, es, ar),
        "ellipse_perimeter_max_diff": max_perim_diff,
    }
