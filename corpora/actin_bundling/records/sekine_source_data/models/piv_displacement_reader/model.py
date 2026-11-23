"""Reads Source-data_Fig5c.xlsx and Source-data_Fig5f.xlsx directly from $BOXGRAPH_SOURCE.

Note: Both spreadsheets label columns 'Vx (nm/s)' and 'Vy (nm/s)', though Fig5c stores unscaled px/s
(~43.25x smaller than plotted nm/s) and Fig5f stores um/s (1000x smaller than plotted nm/s). Outputs
return the exact stored cell values verbatim under the stated column header unit so they compare
directly with the spreadsheets' own AVERAGE/STDEV cells, while Vy/Vx and P-values are scale-invariant.
"""
from __future__ import annotations

import os
import statistics
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

from scipy import stats

_M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"

_SPECS = {
    "fig5c": (
        "Source-data_Fig5c.xlsx",
        "Figure 5c",
        {
            "-40 min": ("B", "C", "D"),
            "-30 min": ("F", "G", "H"),
            "-20 min": ("J", "K", "L"),
            "-10 min": ("N", "O", "P"),
            "0 min": ("R", "S", "T"),
        },
    ),
    "fig5f": (
        "Source-data_Fig5f.xlsx",
        "Figure 5f",
        {
            "Control": ("B", "C", "D"),
            "DAAM RNAi": ("F", "G", "H"),
            "zasp52 RNAi": ("J", "K", "L"),
            "myosin II[DN]": ("N", "O", "P"),
        },
    ),
}


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
    dataset = str(settings.get("dataset", "fig5f"))
    condition = str(settings.get("condition", "Control"))
    fname, sheet, conds = _SPECS[dataset]
    ic, xc, yc = conds[condition]
    rows = _read_sheet(src / fname, sheet)

    vxs, vys = [], []
    for r in sorted(rows):
        if r < 7:
            continue
        if ic in rows[r] and isinstance(rows[r][ic], (int, float)):
            vxs.append(float(rows[r][xc]))
            vys.append(float(rows[r][yc]))

    mx, my = statistics.mean(vxs), statistics.mean(vys)
    return {
        "n_images": float(len(vxs)),
        "vx_mean": mx,
        "vx_sd": statistics.stdev(vxs),
        "vy_mean": my,
        "vy_sd": statistics.stdev(vys),
        "aniso_ratio": my / mx,
        "welch_p": float(stats.ttest_ind(vxs, vys, equal_var=False).pvalue),
        "student_p": float(stats.ttest_ind(vxs, vys, equal_var=True).pvalue),
    }
