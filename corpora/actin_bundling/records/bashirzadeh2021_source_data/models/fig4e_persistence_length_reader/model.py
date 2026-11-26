"""Live reader auditing Fig. 4e and Supp. Fig. 9c bundle persistence length data."""
from __future__ import annotations

import math
import os
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np

_M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def _read_sheet_cells(path: Path, sheet_name: str) -> dict[str, tuple[str | None, str | None]]:
    with zipfile.ZipFile(path) as z:
        sst: list[str] = []
        if "xl/sharedStrings.xml" in z.namelist():
            root = ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in root.findall(_M + "si"):
                sst.append("".join(t.text or "" for t in si.iter(_M + "t")))
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        relmap = {r.get("Id"): r.get("Target") for r in rels}
        rid = next(s.get(_R + "id") for s in wb.iter(_M + "sheet") if s.get("name") == sheet_name)
        target = relmap[rid].lstrip("/")
        if not target.startswith("xl/"):
            target = "xl/" + target
        cells: dict[str, tuple[str | None, str | None]] = {}
        for c in ET.fromstring(z.read(target)).iter(_M + "c"):
            ref = c.get("r")
            t = c.get("t")
            f = c.find(_M + "f")
            v = c.find(_M + "v")
            val = sst[int(v.text)] if (t == "s" and v is not None) else (v.text if v is not None else None)
            ftxt = f.text if f is not None else None
            if val is not None or f is not None:
                cells[ref] = (val, ftxt)
    return cells


def _col_floats(cells: dict[str, tuple[str | None, str | None]], col: str, r1: int, r2: int) -> list[float]:
    out: list[float] = []
    for r in range(r1, r2 + 1):
        entry = cells.get(f"{col}{r}")
        if entry is not None and entry[0] is not None and entry[0].strip() != "":
            out.append(float(entry[0]))
    return out


# Condition index -> (lp_col, len_col, r1, r2, sheet_sem_denom, sd_end_row)
_SPECS = {
    1: ("A", "B", 2, 23, 13, 23),   # [0:1:10] No fascin, 0.1 (A32 = A33/SQRT(13))
    2: ("G", "H", 2, 15, 14, 15),   # [0:3:10] No fascin, 0.3 (G32 = G33/SQRT(14))
    3: ("K", "L", 2, 19, 15, 19),   # [1:1:10] +Fascin 0.1, 0.1 (K32 = K33/SQRT(15), 18 rows)
    4: ("Q", "R", 2, 15, 14, 15),   # [1:2:10] +Fascin 0.1, 0.2 (Q32 = Q33/SQRT(14))
    5: ("X", "Y", 2, 27, 24, 25),   # [1:3:10] +Fascin 0.1, 0.3 (X33 = STDEV(X2:X25), X32 = X33/SQRT(24))
    6: ("AG", None, 3, 17, 15, 17), # [3:3:10] Supp. Fig. 9c (+Fascin 0.3, 0.3; 15 rows AG3:AG17)
}


def run(inputs, settings):
    src = Path(os.environ["BOXGRAPH_SOURCE"])
    f4_path = src / "Fig4e_persistence_length_data.xlsx"
    sd_path = src / "Source_data.xlsx"

    s5 = _read_sheet_cells(f4_path, "Sheet5")
    s5_sd = _read_sheet_cells(sd_path, "Sheet5")
    ppl = _read_sheet_cells(f4_path, "P_PL")
    ppl_sd = _read_sheet_cells(sd_path, "P_PL")
    match_flag = 1.0 if (s5 == s5_sd and ppl == ppl_sd) else 0.0

    cond = int(round(float(inputs["condition"])))
    lp_col, len_col, r1, r2, sheet_denom, sd_r2 = _SPECS[cond]

    vals = _col_floats(s5, lp_col, r1, r2)
    arr = np.array(vals, dtype=float)
    sd_subset = np.array(_col_floats(s5, lp_col, r1, sd_r2), dtype=float)

    uniq = np.array(list(dict.fromkeys(vals)), dtype=float)

    if len_col is not None:
        lengths = np.array(_col_floats(s5, len_col, r1, r2), dtype=float)
        n_with_len = len(lengths)
        n_oor = int(((lengths < 8.0) | (lengths > 20.0)).sum())
    else:
        n_with_len = 0
        n_oor = 0

    sd_full = float(arr.std(ddof=1))
    sd_sheet = float(sd_subset.std(ddof=1))

    return {
        "n_bundles": float(len(arr)),
        "n_unique_bundles": float(len(uniq)),
        "lp_mean": float(arr.mean()),
        "lp_sd": sd_full,
        "lp_sem_true": sd_full / math.sqrt(len(arr)),
        "lp_sem_sheet_formula": sd_sheet / math.sqrt(sheet_denom),
        "lp_median": float(np.median(arr)),
        "lp_dedup_mean": float(uniq.mean()),
        "lp_dedup_sem": float(uniq.std(ddof=1)) / math.sqrt(len(uniq)),
        "n_with_length": float(n_with_len),
        "n_length_out_of_range": float(n_oor),
        "workbook_match_flag": match_flag,
    }
