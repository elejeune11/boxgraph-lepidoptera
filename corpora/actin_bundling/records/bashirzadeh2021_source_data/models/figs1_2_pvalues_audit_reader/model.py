"""Live audit of all 18 Welch t-tests in Figs_1_and_2_Pvalues.xlsx vs Source_data.xlsx."""
from __future__ import annotations

import os
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
from scipy import stats

_M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def _load_workbook_cells(path: Path) -> dict[str, dict[str, tuple[str | None, str | None]]]:
    out: dict[str, dict[str, tuple[str | None, str | None]]] = {}
    with zipfile.ZipFile(path) as z:
        sst: list[str] = []
        if "xl/sharedStrings.xml" in z.namelist():
            root = ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in root.findall(_M + "si"):
                sst.append("".join(t.text or "" for t in si.iter(_M + "t")))
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        relmap = {r.get("Id"): r.get("Target") for r in rels}
        for s in wb.iter(_M + "sheet"):
            name, rid = s.get("name"), s.get(_R + "id")
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
            out[name] = cells
    return out


def _f(cells: dict[str, tuple[str | None, str | None]], ref: str) -> float:
    return float(cells[ref][0])


def _welch_excel(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float, float, float]:
    t_w, p_exact = stats.ttest_ind(x, y, equal_var=False)
    v1, v2 = x.var(ddof=1), y.var(ddof=1)
    n1, n2 = len(x), len(y)
    df_ws = (v1 / n1 + v2 / n2) ** 2 / ((v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1))
    df_ex = float(round(df_ws))
    p_ex = float(2.0 * stats.t.sf(abs(t_w), df_ex))
    return float(t_w), df_ex, float(df_ws), p_ex, float(p_exact)


# Maps test_id (1..18 -> Sheet3..Sheet20) to:
# (pv_sheet_name, pv_rows, pv_col1, pv_col2, sd_sheet, sd_refs1, sd_refs2, aster_row_or_None)
_TESTS = {
    # Sheets 3-5: Fig. 2d Ring (alpha+fascin, 7-12 vs >16) — mistakenly printed on Fig. 1c!
    1:  ("Sheet3",  (1, 2, 3),    "A", "G", "Sheet3", ("N17", "N20", "N23"), ("N19", "N22", "N25"), None),
    2:  ("Sheet4",  (1, 2, 3),    "B", "H", "Sheet3", ("J17", "J20", "J23"), ("J19", "J22", "J25"), None),
    3:  ("Sheet5",  (1, 2, 3),    "C", "I", "Sheet3", ("F17", "F20", "F23"), ("F19", "F22", "F25"), None),
    # Sheets 6-8: Fig. 2e Network (alpha+fascin, 7-12 vs >16) — mistakenly printed on Fig. 1d!
    4:  ("Sheet6",  (6, 7, 8),    "A", "G", "Sheet3", ("P17", "P20", "P23"), ("P19", "P22", "P25"), None),
    5:  ("Sheet7",  (6, 7, 8),    "B", "H", "Sheet3", ("L17", "L20", "L23"), ("L19", "L22", "L25"), None),
    6:  ("Sheet8",  (6, 7, 8),    "C", "I", "Sheet3", ("H17", "H20", "H23"), ("H19", "H22", "H25"), None),
    # Sheets 9-11: Fig. 2f Aster (alpha+fascin, 7-16 vs >16) — mistakenly printed on Fig. 1f!
    7:  ("Sheet9",  (16, 17, 18), "A", "D", "Sheet3", ("M43", "N43", "O43"), ("P43", "Q43", "R43"), (43, 1)),
    8:  ("Sheet10", (16, 17, 18), "B", "E", "Sheet3", ("M44", "N44", "O44"), ("P44", "Q44", "R44"), (44, 2)),
    9:  ("Sheet11", (16, 17, 18), "C", "F", "Sheet3", ("M45", "N45", "O45"), ("P45", "Q45", "R45"), (45, 3)),
    # Sheets 12-14: Fig. 1c Ring (alpha only, 7-12 vs >16) — mistakenly printed on Fig. 2d (with 0.0064 -> 0.0006 typo)!
    10: ("Sheet12", (23, 24, 25), "A", "G", "Sheet2", ("N17", "N20", "N23"), ("N19", "N22", "N25"), None),
    11: ("Sheet13", (23, 24, 25), "B", "H", "Sheet2", ("J17", "J20", "J23"), ("J19", "J22", "J25"), None),
    12: ("Sheet14", (23, 24, 25), "C", "I", "Sheet2", ("F17", "F20", "F23"), ("F19", "F22", "F25"), None),
    # Sheets 15-17: Fig. 1d Network (alpha only, 7-12 vs >16) — mistakenly printed on Fig. 2e!
    13: ("Sheet15", (28, 29, 30), "A", "G", "Sheet2", ("P17", "P20", "P23"), ("P19", "P22", "P25"), None),
    14: ("Sheet16", (28, 29, 30), "B", "H", "Sheet2", ("L17", "L20", "L23"), ("L19", "L22", "L25"), None),
    15: ("Sheet17", (28, 29, 30), "C", "I", "Sheet2", ("H17", "H20", "H23"), ("H19", "H22", "H25"), None),
    # Sheets 18-20: Fig. 1f Peripheral Aster (alpha only, 7-16 vs >16) — mistakenly printed on Fig. 2f!
    16: ("Sheet18", (38, 39, 40), "A", "D", "Sheet2", ("M43", "N43", "O43"), ("P43", "Q43", "R43"), (43, 1)),
    17: ("Sheet19", (38, 39, 40), "B", "E", "Sheet2", ("M44", "N44", "O44"), ("P44", "Q44", "R44"), (44, 2)),
    18: ("Sheet20", (38, 39, 40), "C", "F", "Sheet2", ("M45", "N45", "O45"), ("P45", "Q45", "R45"), (45, 3)),
}


def _sum_range(cells: dict[str, tuple[str | None, str | None]], cols: str, r1: int, r2: int) -> float:
    return sum(_f(cells, f"{c}{r}") for r in range(r1, r2 + 1) for c in cols)


def run(inputs, settings):
    src = Path(os.environ["BOXGRAPH_SOURCE"])
    pv_wb = _load_workbook_cells(src / "Figs_1_and_2_Pvalues.xlsx")
    sd_wb = _load_workbook_cells(src / "Source_data.xlsx")

    tid = int(round(float(inputs["test_id"])))
    sh_name, pv_rows, c1, c2, sd_name, sd_r1, sd_r2, aster_info = _TESTS[tid]

    pv_s1 = pv_wb["Sheet1"]
    sd_sh = sd_wb[sd_name]

    x_pv = np.array([_f(pv_s1, f"{c1}{r}") for r in pv_rows], dtype=float)
    y_pv = np.array([_f(pv_s1, f"{c2}{r}") for r in pv_rows], dtype=float)

    x_sd = np.array([100.0 * _f(sd_sh, ref) for ref in sd_r1], dtype=float)
    y_sd = np.array([100.0 * _f(sd_sh, ref) for ref in sd_r2], dtype=float)

    max_diff = float(max(np.max(np.abs(x_pv - x_sd)), np.max(np.abs(y_pv - y_sd))))

    t_w, df_ex, df_ws, p_ex, p_exact = _welch_excel(x_pv, y_pv)
    _, p_stud = stats.ttest_ind(x_pv, y_pv, equal_var=True)
    p_stored = _f(pv_wb[sh_name], "B12")

    # Recompute with spurious "+1" denominator removed for aster rows (Q43 and O44)
    if aster_info is not None:
        arow, cond = aster_info
        x_fix = x_sd.copy()
        y_fix = y_sd.copy()
        if cond == 1:
            # Q43 = num/(SUM(N3:P3)+1) -> remove +1
            m = re.match(r"^(\d+)/", (sd_sh[f"Q{arow}"][1] or "").strip())
            if m:
                y_fix[1] = 100.0 * float(m.group(1)) / _sum_range(sd_sh, "NOP", 3, 3)
        elif cond == 2:
            # O44 = num/(SUM(J10:L11)+1) -> remove +1
            m = re.match(r"^(\d+)/", (sd_sh[f"O{arow}"][1] or "").strip())
            if m:
                x_fix[2] = 100.0 * float(m.group(1)) / _sum_range(sd_sh, "JKL", 10, 11)
        _, _, _, p_ex_fix, _ = _welch_excel(x_fix, y_fix)
    else:
        p_ex_fix = p_ex

    return {
        "mean_small_pct": float(x_sd.mean()),
        "mean_large_pct": float(y_sd.mean()),
        "t_welch": t_w,
        "df_welch_excel": df_ex,
        "df_welch_exact": df_ws,
        "p_welch_excel_two": p_ex,
        "p_welch_exact_two": p_exact,
        "p_student_two": float(p_stud),
        "p_sheet_stored": p_stored,
        "p_welch_excel_fixed_plus1": p_ex_fix,
        "max_abs_diff_vs_source_data": max_diff,
    }
