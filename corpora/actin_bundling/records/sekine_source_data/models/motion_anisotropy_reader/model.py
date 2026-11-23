"""Reads raw nanocluster tracking coordinates from Source-data_Fig4d.xlsx, Fig5a.xlsx, and Fig5d.xlsx."""
from __future__ import annotations

import math
import os
import statistics
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

from scipy import stats

_M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"

# dataset -> (filename, sheet_name, {condition: (x1, y1, x2, y2, dx, dy)})
_SPECS = {
    "fig4d": (
        "Source-data_Fig4d.xlsx",
        "Figure 4d",
        {
            "1-M": ("C", "D", "F", "G", "I", "J"),
            "2-M": ("N", "O", "Q", "R", "T", "U"),
            "4-M": ("Y", "Z", "AB", "AC", "AE", "AF"),
        },
    ),
    "fig5a": (
        "Source-data_Fig5a.xlsx",
        "Figure 5a",
        {
            "-40 min": ("C", "D", "F", "G", "I", "J"),
            "-30 min": ("N", "O", "Q", "R", "T", "U"),
            "-20 min": ("Y", "Z", "AB", "AC", "AE", "AF"),
            "-10 min": ("AJ", "AK", "AM", "AN", "AP", "AQ"),
            "0 min": ("AU", "AV", "AX", "AY", "BA", "BB"),
        },
    ),
    "fig5d": (
        "Source-data_Fig5d.xlsx",
        "Figure 5d",
        {
            "Control": ("C", "D", "F", "G", "I", "J"),
            "DAAM RNAi": ("N", "O", "Q", "R", "T", "U"),
            "Zasp52 RNAi": ("Y", "Z", "AB", "AC", "AE", "AF"),
            "Myosin II[DN]": ("AJ", "AK", "AM", "AN", "AP", "AQ"),
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
    dataset = str(settings.get("dataset", "fig4d"))
    condition = str(settings.get("condition", "1-M"))
    override_scale = float(settings.get("override_um_per_px", 0.0))

    fname, sheet, conds = _SPECS[dataset]
    x1c, y1c, x2c, y2c, dxc, dyc = conds[condition]
    rows = _read_sheet(src / fname, sheet)

    dxs, dys, dists, implied = [], [], [], []
    for r in sorted(rows):
        if r < 8 or dxc not in rows[r] or dyc not in rows[r]:
            continue
        dx_raw = rows[r][dxc]
        dy_raw = rows[r][dyc]
        if not isinstance(dx_raw, (int, float)) or not isinstance(dy_raw, (int, float)):
            continue
        x1 = float(rows[r][x1c])
        y1 = float(rows[r][y1c])
        x2 = float(rows[r][x2c])
        y2 = float(rows[r][y2c])
        if abs(x2 - x1) > 1e-9:
            implied.append(abs(float(dx_raw) / (x2 - x1)))
        if override_scale > 0.0:
            dx_nm = abs(x2 - x1) * override_scale * 1000.0
            dy_nm = abs(y2 - y1) * override_scale * 1000.0
        else:
            dx_nm = abs(float(dx_raw)) * 1000.0
            dy_nm = abs(float(dy_raw)) * 1000.0
        dxs.append(dx_nm)
        dys.append(dy_nm)
        dists.append(math.hypot(dx_nm, dy_nm))

    mx = statistics.mean(dxs)
    my = statistics.mean(dys)
    return {
        "n_clusters": float(len(dxs)),
        "dx_mean": mx,
        "dx_sd": statistics.stdev(dxs),
        "dy_mean": my,
        "dy_sd": statistics.stdev(dys),
        "aniso_ratio": my / mx,
        "dist_mean": statistics.mean(dists),
        "dist_sd": statistics.stdev(dists),
        "welch_p": float(stats.ttest_ind(dxs, dys, equal_var=False).pvalue),
        "implied_um_per_px": statistics.median(implied) if implied else float("nan"),
    }
