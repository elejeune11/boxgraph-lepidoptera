"""Reads ring thickness vs effective diameter from MOESM14 ('Fig 3e','Fig 3f'), MOESM15 ('Fig 4h'), MOESM20 ('Fig S4d')."""
from __future__ import annotations

import math
import os
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
from scipy.optimize import curve_fit

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


def _pairs(rows, cx: str, cy: str, r0: int = 3) -> list[tuple[float, float]]:
    return [
        (float(rows[r][cx]), float(rows[r][cy]))
        for r in sorted(rows)
        if r >= r0 and cx in rows[r] and cy in rows[r]
    ]


def _pair_set(pts: list[tuple[float, float]]) -> set[tuple[float, float]]:
    return {(round(x, 6), round(y, 6)) for x, y in pts}


def run(inputs, settings):
    src = Path(os.environ["BOXGRAPH_SOURCE"])
    ds = settings.get("dataset", "fig3e_rings")
    d_query = float(inputs["D"])

    r3e = _read_sheet(src / "41567_2022_1924_MOESM14_ESM.xlsx", "Fig 3e")
    r3f = _read_sheet(src / "41567_2022_1924_MOESM14_ESM.xlsx", "Fig 3f")
    pts_3e_discs = _pairs(r3e, "C", "D")
    pts_3e_rings = _pairs(r3e, "G", "H")
    pts_3f_lt = _pairs(r3f, "C", "D")
    pts_3f_gt = _pairs(r3f, "F", "G")

    if ds == "fig3e_discs":
        pts = pts_3e_discs
    elif ds == "fig3e_rings":
        pts = pts_3e_rings
    elif ds == "fig3f_ar_lt11":
        pts = pts_3f_lt
    elif ds == "fig3f_ar_gt11":
        pts = pts_3f_gt
    elif ds == "fig4h_ar_lt11":
        r4h = _read_sheet(src / "41567_2022_1924_MOESM15_ESM.xlsx", "Fig 4h")
        pts = _pairs(r4h, "C", "D")
    elif ds == "edfig4d_all":
        rs4d = _read_sheet(src / "41567_2022_1924_MOESM20_ESM.xlsx", "Fig S4d")
        pts = _pairs(rs4d, "C", "D") + _pairs(rs4d, "F", "G") + _pairs(rs4d, "I", "J")
    else:
        raise ValueError(f"unknown dataset: {ds!r}")

    xs = np.array([p[0] for p in pts], dtype=float)
    ys = np.array([p[1] for p in pts], dtype=float)

    # Log-log OLS power-law fit
    lx, ly = np.log(xs), np.log(ys)
    b_log, log_a = np.polyfit(lx, ly, 1)
    a_log = math.exp(float(log_a))
    y_log_hat = float(log_a) + float(b_log) * lx
    r2_log = float(1.0 - np.sum((ly - y_log_hat) ** 2) / np.sum((ly - np.mean(ly)) ** 2))

    # Nonlinear least squares power-law fit y = a * x^b (reproduces all 5 captions to 4 decimals)
    popt, _ = curve_fit(lambda x, a, b: a * (x ** b), xs, ys, p0=[a_log, float(b_log)])
    a_nl, b_nl = float(popt[0]), float(popt[1])
    y_nl_hat = a_nl * (xs ** b_nl)
    r2_nl = float(1.0 - np.sum((ys - y_nl_hat) ** 2) / np.sum((ys - np.mean(ys)) ** 2))

    pset = _pair_set(pts)
    s3e = _pair_set(pts_3e_discs + pts_3e_rings)
    s3f = _pair_set(pts_3f_lt)

    return {
        "thickness_nls": a_nl * (d_query ** b_nl),
        "a_nls": a_nl,
        "b_nls": b_nl,
        "r2_nls": r2_nl,
        "a_logols": a_log,
        "b_logols": float(b_log),
        "r2_logols": r2_log,
        "n_droplets": float(len(pts)),
        "diam_min": float(xs.min()),
        "diam_max": float(xs.max()),
        "n_shared_fig3e": float(len(pset & s3e)),
        "n_shared_fig3f_lt11": float(len(pset & s3f)),
    }
