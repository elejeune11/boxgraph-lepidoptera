"""Reads MOESM16 ('Fig S5a', 'Fig S5b') and MOESM21 ('Fig S8a') across actin:VASP ratios 1:20, 1:10, 1:7.5."""
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


def _col_vals(rows, col: str, r0: int) -> np.ndarray:
    return np.array([float(rows[r][col]) for r in sorted(rows) if r >= r0 and col in rows[r]], dtype=float)


def _ratio_key(ratio: float) -> str:
    options = [("1:20", 1.0 / 20.0), ("1:10", 1.0 / 10.0), ("1:7.5", 1.0 / 7.5)]
    return min(options, key=lambda kv: abs(ratio - kv[1]))[0]


def run(inputs, settings):
    src = Path(os.environ["BOXGRAPH_SOURCE"])
    rk = _ratio_key(float(inputs["ratio"]))

    r5a = _read_sheet(src / "41567_2022_1924_MOESM16_ESM.xlsx", "Fig S5a")
    r5b = _read_sheet(src / "41567_2022_1924_MOESM16_ESM.xlsx", "Fig S5b")
    rs8 = _read_sheet(src / "41567_2022_1924_MOESM21_ESM.xlsx", "Fig S8a")

    arr_5a = {
        "1:20": _col_vals(r5a, "C", 4),
        "1:10": _col_vals(r5a, "D", 4),
        "1:7.5": _col_vals(r5a, "E", 4),
    }
    next_map = {"1:20": "1:10", "1:10": "1:7.5", "1:7.5": "1:20"}
    a = arr_5a[rk]
    _, p_stud = stats.ttest_ind(a, arr_5a[next_map[rk]], equal_var=True)
    q25, med, q75 = np.percentile(a, [25.0, 50.0, 75.0])

    int_5b = {
        "1:20": _col_vals(r5b, "C", 3),
        "1:10": _col_vals(r5b, "F", 3),
        "1:7.5": _col_vals(r5b, "I", 3),
    }
    s8_rep2 = {
        "1:20": _col_vals(rs8, "C", 5),
        "1:10": _col_vals(rs8, "E", 5),
        "1:7.5": _col_vals(rs8, "G", 5),
    }
    s8_rep3 = {
        "1:20": _col_vals(rs8, "J", 5),
        "1:10": _col_vals(rs8, "L", 5),
        "1:7.5": _col_vals(rs8, "N", 5),
    }

    return {
        "n_droplets_5a": float(len(a)),
        "ar_median": float(med),
        "ar_mean": float(a.mean()),
        "ar_sd": float(a.std(ddof=1)),
        "ar_q25": float(q25),
        "ar_q75": float(q75),
        "ar_max": float(a.max()),
        "p_val_student_vs_next": float(p_stud),
        "n_droplets_5b": float(len(int_5b[rk])),
        "total_5b": float(sum(len(v) for v in int_5b.values())),
        "intensity_mean_5b": float(int_5b[rk].mean()),
        "n_droplets_s8a_rep2": float(len(s8_rep2[rk])),
        "total_s8a_rep2": float(sum(len(v) for v in s8_rep2.values())),
        "n_droplets_s8a_rep3": float(len(s8_rep3[rk])),
        "total_s8a_rep3": float(sum(len(v) for v in s8_rep3.values())),
    }
