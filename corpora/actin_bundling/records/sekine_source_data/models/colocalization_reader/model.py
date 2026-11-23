"""Reads Source-data_Fig3d.xlsx (Stage 15) and Source-data_SupplementaryFig5d.xlsx (Stage 16) from $BOXGRAPH_SOURCE."""
from __future__ import annotations

import os
import statistics
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

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


def _load_groups(src: Path, stage: int) -> dict[str, list[float]]:
    if stage == 15:
        rows = _read_sheet(src / "Source-data_Fig3d.xlsx", "Figure 3d")
    elif stage == 16:
        rows = _read_sheet(src / "Source-data_SupplementaryFig5d.xlsx", "Supplementary Figure 5d")
    else:
        raise ValueError(f"unsupported stage: {stage}")
    groups: dict[str, list[float]] = {}
    cur = None
    for r in sorted(rows):
        if r < 7:
            continue
        if "B" in rows[r] and isinstance(rows[r]["B"], str):
            cur = rows[r]["B"].strip()
        if "D" in rows[r] and isinstance(rows[r]["D"], (int, float)) and cur is not None:
            groups.setdefault(cur, []).append(float(rows[r]["D"]))
    return groups


def run(inputs, settings):
    src = Path(os.environ["BOXGRAPH_SOURCE"])
    stage = int(settings.get("stage", 16))
    protein = str(settings.get("protein", "RhoGEF2::mCherry"))
    control_stage = int(settings.get("control_stage", stage))

    groups = _load_groups(src, stage)
    ctrl_groups = groups if control_stage == stage else _load_groups(src, control_stage)
    vals = groups[protein]
    gfp = ctrl_groups["GFP"]

    if protein == "GFP" and control_stage == stage:
        pw, pt = 1.0, 1.0
    else:
        pw = float(stats.ttest_ind(vals, gfp, equal_var=False).pvalue)
        pt = float(stats.ttest_ind(vals, gfp, equal_var=True).pvalue)

    return {
        "n_rois": float(len(vals)),
        "r_mean": statistics.mean(vals),
        "r_sd": statistics.stdev(vals),
        "r_median": statistics.median(vals),
        "welch_p_vs_gfp": pw,
        "student_p_vs_gfp": pt,
    }
