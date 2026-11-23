"""Reads raw per-tube fusion events in Source-data_Fig2b.xlsx directly from $BOXGRAPH_SOURCE."""
from __future__ import annotations

import os
import statistics
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

from scipy import stats

_M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def _col2num(col: str) -> int:
    n = 0
    for ch in col:
        n = n * 26 + ord(ch) - 64
    return n


def _num2col(n: int) -> str:
    s = ""
    while n > 0:
        n, r = divmod(n - 1, 26)
        s = chr(65 + r) + s
    return s


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


def run(inputs, settings):
    src = Path(os.environ["BOXGRAPH_SOURCE"])
    t_target = int(settings.get("time_min", 0))
    rows = _read_sheet(src / "Source-data_Fig2b.xlsx", "Figure 2b")

    # Raw Angle in cols I-EJ is measured from circumferential axis (0 deg = circumferential, 90 deg = longitudinal),
    # so theta relative to tube longitudinal axis is 90 - min(Angle, 180 - Angle).
    circum, diag, longit = [], [], []
    for col, hdr in rows[6].items():
        if _col2num(col) < _col2num("I") or not isinstance(hdr, str) or "(" not in hdr:
            continue
        t_block = int(hdr.split("(")[1].split(")")[0].replace("min", "").strip())
        if t_block != t_target:
            continue
        cidx = _col2num(col)
        dur_col = _num2col(cidx + 3)
        ang_col = _num2col(cidx + 4)
        for r in sorted(rows):
            if r < 8:
                continue
            if dur_col in rows[r] and ang_col in rows[r]:
                dur = rows[r][dur_col]
                ang = rows[r][ang_col]
                if isinstance(dur, (int, float)) and isinstance(ang, (int, float)):
                    folded = float(ang) if float(ang) <= 90.0 else 180.0 - float(ang)
                    theta = 90.0 - folded
                    if theta > 60.0:
                        circum.append(float(dur))
                    elif theta >= 30.0:
                        diag.append(float(dur))
                    else:
                        longit.append(float(dur))

    kw_p = float(stats.kruskal(circum, diag, longit).pvalue)
    mwu_p = float(stats.mannwhitneyu(circum, diag, alternative="two-sided").pvalue)
    return {
        "n_circum": float(len(circum)),
        "n_diag": float(len(diag)),
        "n_long": float(len(longit)),
        "mean_dur_circum": statistics.mean(circum),
        "mean_dur_diag": statistics.mean(diag),
        "mean_dur_long": statistics.mean(longit),
        "median_dur_circum": statistics.median(circum),
        "median_dur_diag": statistics.median(diag),
        "median_dur_long": statistics.median(longit),
        "kw_p": kw_p,
        "wilcox_circum_vs_diag_unadj_p": mwu_p,
        "wilcox_circum_vs_diag_bonf_p": min(1.0, 3.0 * mwu_p),
    }
