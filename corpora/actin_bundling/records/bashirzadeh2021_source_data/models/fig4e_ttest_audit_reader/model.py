"""Live statistical audit of Fig. 4e and Supp. Fig. 9c persistence length comparisons."""
from __future__ import annotations

import os
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
from scipy import stats

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


def _col(cells: dict[str, float], col: str, r1: int, r2: int) -> list[float]:
    return [cells[f"{col}{r}"] for r in range(r1, r2 + 1) if f"{col}{r}" in cells]


def run(inputs, settings):
    src = Path(os.environ["BOXGRAPH_SOURCE"])
    cells = _read_sheet_floats(src / "Fig4e_persistence_length_data.xlsx", "Sheet5")

    groups = {
        "0_1_10": _col(cells, "A", 2, 23),
        "0_3_10": _col(cells, "G", 2, 15),
        "1_1_10": _col(cells, "K", 2, 19),
        "1_3_10": _col(cells, "X", 2, 27),
        "3_3_10": _col(cells, "AG", 3, 17),
    }

    # Pair ordering matches P_PL (Variable 1 minus Variable 2):
    # 1: [1:1:10] vs [0:1:10]  (0-0.1 vs 0.1-0.1)
    # 2: [1:3:10] vs [1:1:10]  (0.1-0.1 vs 0.1-0.3)
    # 3: [1:3:10] vs [0:3:10]  (0-0.3 vs 0.1-0.3)
    # 4: [3:3:10] vs [1:3:10]  (Supp. Fig. 9c: fascin/actin 0.3 vs 0.1 at alpha=0.3)
    pair_map = {
        1: ("1_1_10", "0_1_10"),
        2: ("1_3_10", "1_1_10"),
        3: ("1_3_10", "0_3_10"),
        4: ("3_3_10", "1_3_10"),
    }
    k1, k2 = pair_map[int(round(float(inputs["pair"])))]
    v1_list, v2_list = groups[k1], groups[k2]
    a = np.array(v1_list, dtype=float)
    b = np.array(v2_list, dtype=float)

    t_w, p_w_exact = stats.ttest_ind(a, b, equal_var=False)
    s1, s2 = a.var(ddof=1), b.var(ddof=1)
    n1, n2 = len(a), len(b)
    df_ws = (s1 / n1 + s2 / n2) ** 2 / ((s1 / n1) ** 2 / (n1 - 1) + (s2 / n2) ** 2 / (n2 - 1))
    # Excel Analysis ToolPak rounds Welch-Satterthwaite df to nearest integer.
    df_excel = float(round(df_ws))
    p_excel_one = float(stats.t.sf(abs(t_w), df_excel))
    p_excel_two = 2.0 * p_excel_one

    _, p_stud = stats.ttest_ind(a, b, equal_var=True)
    _, p_mw = stats.mannwhitneyu(a, b, alternative="two-sided")
    _, p_log = stats.ttest_ind(np.log(a), np.log(b), equal_var=False)

    a_u = np.array(list(dict.fromkeys(v1_list)), dtype=float)
    b_u = np.array(list(dict.fromkeys(v2_list)), dtype=float)
    _, p_dedup = stats.ttest_ind(a_u, b_u, equal_var=False)

    return {
        "t_welch": float(t_w),
        "df_welch_excel": df_excel,
        "df_welch_exact": float(df_ws),
        "p_welch_excel_two": p_excel_two,
        "p_welch_excel_one": p_excel_one,
        "p_welch_exact_two": float(p_w_exact),
        "p_student_two": float(p_stud),
        "p_mann_whitney_two": float(p_mw),
        "p_log_welch_two": float(p_log),
        "p_dedup_welch_two": float(p_dedup),
    }
