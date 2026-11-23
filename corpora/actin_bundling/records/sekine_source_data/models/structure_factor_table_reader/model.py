"""Reads Source-data_SupplementaryFig7c.xlsx (and checks equivalence with stray sheet in Fig5f.xlsx)."""
from __future__ import annotations

import math
import os
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

_M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"

_ROW_BY_COND = {"1-M": 9, "3-M": 10, "4-M": 11, "a4-1": 12}


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


def run(inputs, settings):
    src = Path(os.environ["BOXGRAPH_SOURCE"])
    cond = str(settings.get("condition", "1-M"))
    r7c = _read_sheet(src / "Source-data_SupplementaryFig7c.xlsx", "Supplementary Figure 7c")
    r5f_dup = _read_sheet(src / "Source-data_Fig5f.xlsx", "Supplementary Figure 6e")
    rn = _ROW_BY_COND[cond]
    if r7c[rn] != r5f_dup[rn]:
        raise ValueError(f"stray Fig5f sheet differs from SuppFig7c at row {rn}")
    A = float(r7c[rn]["C"])
    mu_2pi = float(r7c[rn]["D"])
    sigma = float(r7c[rn]["E"])
    interval_stored = float(r7c[rn]["F"])
    return {
        "A": A,
        "mu_over_2pi": mu_2pi,
        "sigma": sigma,
        "interval_stored": interval_stored,
        "interval_from_mu": 1000.0 / mu_2pi,
        "peak_height": A / (math.sqrt(2.0 * math.pi) * sigma),
    }
