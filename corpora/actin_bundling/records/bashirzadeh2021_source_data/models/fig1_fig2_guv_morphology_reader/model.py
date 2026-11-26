"""Live reader auditing Source_data.xlsx Sheet2 (Fig. 1) and Sheet3 (Fig. 2)."""
from __future__ import annotations

import math
import os
import re
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


def _num(cells: dict[str, tuple[str | None, str | None]], ref: str) -> float:
    entry = cells.get(ref)
    if entry is None or entry[0] is None:
        return 0.0
    return float(entry[0])


def _row_triple(cells: dict[str, tuple[str | None, str | None]], cols: str, row: int) -> tuple[float, float, float]:
    return (_num(cells, f"{cols[0]}{row}"), _num(cells, f"{cols[1]}{row}"), _num(cells, f"{cols[2]}{row}"))


# Per-condition 3 replicates: each entry is (cols, (row_small, row_med, row_large))
_REP_BLOCKS = {
    1: [("NOP", (1, 2, 3)), ("NOP", (5, 6, 7)), ("QRS", (11, 12, 13))],  # 0.1
    2: [("JKL", (4, 5, 6)), ("JKL", (7, 8, 9)), ("JKL", (10, 11, 12))],  # 0.2
    3: [("FGH", (4, 5, 6)), ("FGH", (7, 8, 9)), ("FGH", (10, 11, 12))],  # 0.3
}


def run(inputs, settings):
    src = Path(os.environ["BOXGRAPH_SOURCE"])
    dataset = settings.get("dataset", "fig1_sheet2")
    sheet = "Sheet2" if dataset == "fig1_sheet2" else "Sheet3"
    cells = _read_sheet_cells(src / "Source_data.xlsx", sheet)

    cond = int(round(float(inputs["condition"])))
    sidx = int(round(float(inputs["size_bin"]))) - 1  # 0=7-12, 1=12-16, 2=>16

    reps = _REP_BLOCKS[cond]
    rep_totals: list[float] = []
    size_total = 0.0
    s7to16_total = 0.0
    ring_pcts, rn_pcts, net_pcts = [], [], []

    for cols, rows in reps:
        t0 = _row_triple(cells, cols, rows[0])
        t1 = _row_triple(cells, cols, rows[1])
        t2 = _row_triple(cells, cols, rows[2])
        rep_totals.append(sum(t0) + sum(t1) + sum(t2))
        s7to16_total += sum(t0) + sum(t1)
        chosen = (t0, t1, t2)[sidx]
        tot = sum(chosen)
        size_total += tot
        ring_pcts.append(100.0 * chosen[0] / tot)
        rn_pcts.append(100.0 * chosen[1] / tot)
        net_pcts.append(100.0 * chosen[2] / tot)

    def _mean_sem(xs: list[float]) -> tuple[float, float]:
        a = np.array(xs, dtype=float)
        return float(a.mean()), float(a.std(ddof=1) / math.sqrt(len(a)))

    ring_m, ring_se = _mean_sem(ring_pcts)
    rn_m, rn_se = _mean_sem(rn_pcts)
    net_m, net_se = _mean_sem(net_pcts)

    # Pooled cumulative probabilities across ALL 3 concentrations for this size_bin
    c_ring, c_rn, c_net = 0.0, 0.0, 0.0
    for c_k in (1, 2, 3):
        for cols, rows in _REP_BLOCKS[c_k]:
            r_v, rn_v, n_v = _row_triple(cells, cols, rows[sidx])
            c_ring += r_v
            c_rn += rn_v
            c_net += n_v
    c_tot = c_ring + c_rn + c_net

    # Spreadsheet A8:A10 cached cumulative ring probability (omits J11 for medium bin)
    a8_ring = _num(cells, f"A{8 + sidx}")

    # Aster percentages (Fig. 1f / Fig. 2f): rows 43 (0.1), 44 (0.2), 45 (0.3)
    arow = 42 + cond
    a7to16_sheet = [100.0 * _num(cells, f"{col}{arow}") for col in ("M", "N", "O")]
    agt16_sheet = [100.0 * _num(cells, f"{col}{arow}") for col in ("P", "Q", "R")]
    a7to16_m, a7to16_se = _mean_sem(a7to16_sheet)
    agt16_m, agt16_se = _mean_sem(agt16_sheet)

    # Parse numerators from P43:R45 formulas to audit true count sum and +1 denominator fix
    gt16_nums: list[float] = []
    # Note: in M43:R45 for cond=1 (0.1), rep ordering in columns P,Q,R is rep2(N7), rep1(N3), rep3(Q13)
    gt16_denoms_true = (
        [
            sum(_row_triple(cells, "NOP", 7)),
            sum(_row_triple(cells, "NOP", 3)),
            sum(_row_triple(cells, "QRS", 13)),
        ]
        if cond == 1
        else [sum(_row_triple(cells, cols, rows[2])) for cols, rows in reps]
    )
    agt16_fixed = []
    for col, denom in zip(("P", "Q", "R"), gt16_denoms_true):
        ftxt = cells[f"{col}{arow}"][1] or ""
        m = re.match(r"^(\d+)/", ftxt.strip())
        num = float(m.group(1)) if m else 0.0
        gt16_nums.append(num)
        agt16_fixed.append(100.0 * num / denom)

    # Fig. 1g (rows 48, 49, 50 for cond 1, 2, 3): M=Lumen, N=Membrane (Peripheral)
    grow = 47 + cond
    lumen_pct = _num(cells, f"M{grow}")
    periph_pct = _num(cells, f"N{grow}")
    g_formula = cells[f"M{grow}"][1] or ""
    mg = re.match(r"^\d+/(\d+)\*100$", g_formula.strip())
    g_denom = float(mg.group(1)) if mg else 0.0

    return {
        "n_guv_size_total": float(size_total),
        "n_guv_7to16_total": float(s7to16_total),
        "n_guv_rep1_total": float(rep_totals[0]),
        "n_guv_rep2_total": float(rep_totals[1]),
        "n_guv_rep3_total": float(rep_totals[2]),
        "ring_pct_mean": ring_m,
        "ring_pct_sem": ring_se,
        "ring_net_pct_mean": rn_m,
        "ring_net_pct_sem": rn_se,
        "net_pct_mean": net_m,
        "net_pct_sem": net_se,
        "cumul_ring_prob_correct": float(c_ring / c_tot),
        "cumul_ring_prob_sheet_a8": float(a8_ring),
        "cumul_rn_prob_correct": float(c_rn / c_tot),
        "cumul_net_prob_correct": float(c_net / c_tot),
        "aster_7to16_pct_mean": a7to16_m,
        "aster_7to16_pct_sem": a7to16_se,
        "aster_gt16_pct_mean": agt16_m,
        "aster_gt16_pct_sem": agt16_se,
        "aster_gt16_pct_mean_fixed": float(np.mean(agt16_fixed)),
        "aster_gt16_count_sum": float(sum(gt16_nums)),
        "fig1g_lumen_pct": float(lumen_pct),
        "fig1g_periph_pct": float(periph_pct),
        "fig1g_denom_count": g_denom,
    }
