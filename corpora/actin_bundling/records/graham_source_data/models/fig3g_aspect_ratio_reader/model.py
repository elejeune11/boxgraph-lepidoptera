"""Reads droplet aspect ratios for shells, rings, and discs from MOESM14 ('Fig 3g')."""
from __future__ import annotations

import os
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
from scipy import stats

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
    rows = _read_sheet(src / "41567_2022_1924_MOESM14_ESM.xlsx", "Fig 3g")
    cat = int(round(float(inputs["category"])))
    cols = {1: "C", 2: "D", 3: "E"}
    arrays = {
        k: np.array([float(rows[r][col]) for r in sorted(rows) if r >= 3 and col in rows[r]], dtype=float)
        for k, col in cols.items()
    }
    arr = arrays[cat]
    # Compare 1 vs 2 (shell vs ring), 2 vs 3 (ring vs disc), 3 vs 1 (disc vs shell)
    other = arrays[{1: 2, 2: 3, 3: 1}[cat]]
    _, p_stud = stats.ttest_ind(arr, other, equal_var=True)
    _, p_welch = stats.ttest_ind(arr, other, equal_var=False)
    q25, med, q75 = np.percentile(arr, [25.0, 50.0, 75.0])
    return {
        "n_droplets": float(len(arr)),
        "ar_median": float(med),
        "ar_mean": float(arr.mean()),
        "ar_sd": float(arr.std(ddof=1)),
        "ar_q25": float(q25),
        "ar_q75": float(q75),
        "ar_max": float(arr.max()),
        "p_val_student_vs_next": float(p_stud),
        "p_val_welch_vs_next": float(p_welch),
    }
