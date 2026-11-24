"""Reads PEG-dependent droplet aspect ratios from MOESM22 ('Fig S9e')."""
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
    rows = _read_sheet(src / "41567_2022_1924_MOESM22_ESM.xlsx", "Fig S9e")
    peg_in = float(inputs["peg"])
    peg_key = min([3, 5, 10, 15], key=lambda p: abs(peg_in - p))
    col_map = {3: "B", 5: "C", 10: "D", 15: "E"}
    arrays = {
        p: np.array([float(rows[r][c]) for r in sorted(rows) if r >= 3 and c in rows[r]], dtype=float)
        for p, c in col_map.items()
    }
    arr = arrays[peg_key]
    other = arrays[{3: 5, 5: 10, 10: 15, 15: 3}[peg_key]]
    _, p_stud = stats.ttest_ind(arr, other, equal_var=True)
    _, p_welch = stats.ttest_ind(arr, other, equal_var=False)
    q25, med, q75 = np.percentile(arr, [25.0, 50.0, 75.0])
    return {
        "n_droplets": float(len(arr)),
        "ar_mean": float(arr.mean()),
        "ar_sd": float(arr.std(ddof=1)),
        "ar_median": float(med),
        "ar_q25": float(q25),
        "ar_q75": float(q75),
        "ar_max": float(arr.max()),
        "n_gt_15": float((arr > 15.0).sum()),
        "p_val_student_vs_next": float(p_stud),
        "p_val_welch_vs_next": float(p_welch),
    }
