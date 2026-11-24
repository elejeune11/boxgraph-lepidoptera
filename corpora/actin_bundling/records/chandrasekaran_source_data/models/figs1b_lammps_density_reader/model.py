"""Reads sheet 'FigureS1B' (LAMMPS radial VASP and crowder number density profiles)."""
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
def _load_s1b(src_str: str):
    rows = _read_sheet(Path(src_str) / _XLSX, "FigureS1B")
    rs, vd, cd = [], [], []
    for rn in range(2, 86):
        rs.append(_f(rows[rn]["A"]))
        vd.append(_f(rows[rn]["B"]))
        cd.append(_f(rows[rn]["C"]))
    r_iface = rs[-1]
    for i in range(len(rs) - 1):
        if vd[i] >= 0.1 >= vd[i + 1]:
            r_iface = rs[i] + (0.1 - vd[i]) / (vd[i + 1] - vd[i]) * (rs[i + 1] - rs[i])
            break
    v_plat = statistics.mean(v for r, v in zip(rs, vd) if 6.0 <= r <= 14.0)
    c_plat = statistics.mean(c for r, c in zip(rs, cd) if 22.0 <= r <= 30.0)
    return rs, vd, cd, r_iface, v_plat, c_plat


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
    rs, vd, cd, r_iface, v_plat, c_plat = _load_s1b(os.environ["BOXGRAPH_SOURCE"])
    rq = float(inputs["r"])
    return {
        "vasp_density": _interp(rs, vd, rq),
        "crowder_density": _interp(rs, cd, rq),
        "interface_radius_0p1": r_iface,
        "vasp_interior_plateau": v_plat,
        "crowder_bulk_plateau": c_plat,
    }
