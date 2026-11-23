"""Reads Source-data_Fig1ef.xlsx and Source-data_SupplementaryFig1bc.xlsx directly from $BOXGRAPH_SOURCE."""
from __future__ import annotations

import math
import os
import statistics
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

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


def _col_floats(rows: dict[int, dict[str, object]], col: str, start_row: int = 7) -> list[float]:
    return [
        float(rows[r][col])
        for r in sorted(rows)
        if r >= start_row and col in rows[r] and isinstance(rows[r][col], (int, float))
    ]


def _ks_2samp_r_asymp(x: list[float], y: list[float]) -> float:
    """Two-sided two-sample KS statistic and R's asymptotic Kolmogorov P-value."""
    sx = sorted(x)
    sy = sorted(y)
    nx, ny = len(sx), len(sy)
    i = j = 0
    d_max = 0.0
    while i < nx and j < ny:
        if sx[i] < sy[j]:
            v = sx[i]
            while i < nx and sx[i] == v:
                i += 1
        elif sy[j] < sx[i]:
            v = sy[j]
            while j < ny and sy[j] == v:
                j += 1
        else:
            v = sx[i]
            while i < nx and sx[i] == v:
                i += 1
            while j < ny and sy[j] == v:
                j += 1
        d = abs(i / nx - j / ny)
        if d > d_max:
            d_max = d
    z = d_max * math.sqrt(nx * ny / (nx + ny))
    if z == 0.0:
        return 1.0
    p = 2.0 * sum(((-1) ** (k - 1)) * math.exp(-2.0 * (k * k) * (z * z)) for k in range(1, 100))
    return max(0.0, min(1.0, p))


def run(inputs, settings):
    src = Path(os.environ["BOXGRAPH_SOURCE"])
    stage = int(settings.get("stage", 15))
    r1ef = _read_sheet(src / "Source-data_Fig1ef.xlsx", "Figure 1e, f")
    rs1bc = _read_sheet(src / "Source-data_SupplementaryFig1bc.xlsx", "Supplementary Figures 1b, c")

    ar15 = _col_floats(r1ef, "I")
    ar16 = _col_floats(r1ef, "U")
    maj15 = _col_floats(r1ef, "E")
    maj16 = _col_floats(r1ef, "Q")
    min15 = _col_floats(r1ef, "F")
    min16 = _col_floats(r1ef, "R")

    if stage == 15:
        area = _col_floats(r1ef, "D")
        major = maj15
        minor = min15
        angle = _col_floats(r1ef, "G")
        ar = ar15
        major_dup = _col_floats(rs1bc, "E")
    elif stage == 16:
        area = _col_floats(r1ef, "P")
        major = maj16
        minor = min16
        angle = _col_floats(r1ef, "S")
        ar = ar16
        major_dup = _col_floats(rs1bc, "Q")
    else:
        raise ValueError(f"unsupported stage: {stage}")

    folded = [a if a <= 90.0 else 180.0 - a for a in angle]
    dup_diff = max(abs(a - b) for a, b in zip(major, major_dup)) if len(major) == len(major_dup) else float("nan")

    return {
        "n_clusters": float(len(major)),
        "area_mean": statistics.mean(area),
        "area_sd": statistics.stdev(area),
        "major_mean": statistics.mean(major),
        "major_sd": statistics.stdev(major),
        "minor_mean": statistics.mean(minor),
        "minor_sd": statistics.stdev(minor),
        "ar_mean": statistics.mean(ar),
        "ar_sd": statistics.stdev(ar),
        "angle_median_folded": statistics.median(folded),
        "ks_p_ar": _ks_2samp_r_asymp(ar15, ar16),
        "ks_p_major": _ks_2samp_r_asymp(maj15, maj16),
        "ks_p_minor": _ks_2samp_r_asymp(min15, min16),
        "dup_s1bc_max_diff": dup_diff,
    }
