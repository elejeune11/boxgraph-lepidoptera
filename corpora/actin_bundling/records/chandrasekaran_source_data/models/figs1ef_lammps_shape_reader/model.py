"""Reads sheet 'FigureS1EF' (LAMMPS droplet asphericity A and aspect ratio c/a at 6 and 20 filaments)."""
from __future__ import annotations

import functools
import os
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
def _load_s1ef(src_str: str):
    rows = _read_sheet(Path(src_str) / _XLSX, "FigureS1EF")
    return {
        _f(rows[rn]["A"]): (_f(rows[rn]["B"]), _f(rows[rn]["C"]))
        for rn in (2, 3)
    }


def run(inputs, settings):
    data = _load_s1ef(os.environ["BOXGRAPH_SOURCE"])
    nq = float(inputs["n_fil"])
    key = min(data, key=lambda k: abs(k - nq))
    asph, ar = data[key]
    # Since A = (c^2 / (a*b))^(1/3) and ar = c/a, b/a = (c/a)^2 / A^3
    b_over_a = (ar * ar) / (asph ** 3)
    return {
        "asphericity": asph,
        "aspect_ratio": ar,
        "implied_b_over_a": b_over_a,
    }
