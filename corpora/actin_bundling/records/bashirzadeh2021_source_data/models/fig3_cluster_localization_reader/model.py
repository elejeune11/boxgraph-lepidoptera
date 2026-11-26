"""Live reader auditing Source_data.xlsx Sheet4 (Figs. 3c and 3d cluster localization)."""
from __future__ import annotations

import os
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

_M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def _read_sheet_floats(path: Path, sheet_name: str) -> dict[str, float]:
    with zipfile.ZipFile(path) as z:
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        relmap = {r.get("Id"): r.get("Target") for r in rels}
        rid = next(s.get(_R + "id") for s in wb.iter(_M + "sheet") if s.get("name") == sheet_name)
        target = relmap[rid].lstrip("/")
        if not target.startswith("xl/"):
            target = "xl/" + target
        out: dict[str, float] = {}
        for c in ET.fromstring(z.read(target)).iter(_M + "c"):
            ref = c.get("r")
            t = c.get("t")
            v = c.find(_M + "v")
            if t is None and v is not None and v.text is not None:
                out[ref] = float(v.text)
    return out


def run(inputs, settings):
    src = Path(os.environ["BOXGRAPH_SOURCE"])
    c = _read_sheet_floats(src / "Source_data.xlsx", "Sheet4")

    cat = int(round(float(inputs["category"])))

    # Summary counts plotted in Figs. 3c and 3d (Sheet4 J4:L4, J6:L6, N5:P5)
    f3d_16_24 = (c["J4"], c["K4"], c["L4"])      # 11, 34, 16 (n=61)
    f3d_gt24 = (c["J6"], c["K6"], c["L6"])       # 15, 7, 4 (n=26)
    f3c_fascin = tuple(a + b for a, b in zip(f3d_16_24, f3d_gt24))  # 26, 41, 20 (n=87)
    f3c_nofascin = (c["N5"], c["O5"], c["P5"])   # 3, 50, 27 (n=80)

    # Stale raw table sums in B4:G6 (16_25 and >25)
    raw_16_25 = (
        c["B4"] + c["B5"] + c["B6"],  # 11
        c["D4"] + c["D5"] + c["D6"],  # 42 (vs K4=34!)
        c["F4"] + c["F5"] + c["F6"],  # 9  (vs L4=16!)
    )
    raw_gt25 = (
        c["C4"] + c["C5"] + c["C6"],  # 6
        c["E4"] + c["E5"] + c["E6"],  # 3
        c["G4"] + c["G5"] + c["G6"],  # 2
    )

    table = {
        1: (f3d_16_24, sum(f3d_16_24) + sum(f3d_gt24)),
        2: (f3d_gt24, sum(f3d_16_24) + sum(f3d_gt24)),
        3: (f3c_fascin, sum(f3c_fascin) + sum(f3c_nofascin)),
        4: (f3c_nofascin, sum(f3c_fascin) + sum(f3c_nofascin)),
        5: (raw_16_25, sum(raw_16_25) + sum(raw_gt25)),
        6: (raw_gt25, sum(raw_16_25) + sum(raw_gt25)),
    }
    (n_center, n_pole, n_none), n_comb = table[cat]
    tot = n_center + n_pole + n_none

    return {
        "n_guvs": float(tot),
        "n_guvs_combined": float(n_comb),
        "center_pct": float(100.0 * n_center / tot),
        "periphery_pct": float(100.0 * n_pole / tot),
        "no_cluster_pct": float(100.0 * n_none / tot),
    }
