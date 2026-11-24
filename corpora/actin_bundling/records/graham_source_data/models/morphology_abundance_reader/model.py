"""Reads morphology count tables from MOESM13 ('Fig 2d'), MOESM14 ('Fig 3d'), MOESM17 ('Fig 6d'), and MOESM20 ('Fig S4c')."""
from __future__ import annotations

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
        rows: dict[int, dict[str, object]] = {}
        for c in ET.fromstring(z.read(target)).iter(_M + "c"):
            ref = c.get("r")
            col = "".join(ch for ch in ref if ch.isalpha())
            rn = int("".join(ch for ch in ref if ch.isdigit()))
            kind, v = c.get("t"), c.find(_M + "v")
            if v is not None and kind is None:
                rows.setdefault(rn, {})[col] = float(v.text)
    return rows


def run(inputs, settings):
    src = Path(os.environ["BOXGRAPH_SOURCE"])
    ds = settings.get("dataset", "fig2d_1to20")
    method = settings.get("method", "pooled")
    cat = int(round(float(inputs["category"])))

    if ds.startswith("fig2d_"):
        rows = _read_sheet(src / "41567_2022_1924_MOESM13_ESM.xlsx", "Fig 2d")
        col_map = {
            "fig2d_1to20": ["N", "O", "P"],
            "fig2d_1to10": ["D", "E", "F"],
            "fig2d_1to7p5": ["I", "J", "K"],
        }
        cols = col_map[ds]
        row_for_cat = {1: 6, 2: 7, 3: 8, 4: 9}[cat]
        rep_tots = [sum(float(rows[r][c]) for r in (6, 7, 8, 9)) for c in cols]
        pcts = [100.0 * float(rows[row_for_cat][c]) / tot for c, tot in zip(cols, rep_tots)]

    elif ds.startswith("fig3d_"):
        rows = _read_sheet(src / "41567_2022_1924_MOESM14_ESM.xlsx", "Fig 3d")
        row_ranges = {"fig3d_1to20": (7, 14), "fig3d_1to10": (20, 27), "fig3d_1to7p5": (32, 37)}
        r0, r1 = row_ranges[ds]
        blocks = [("C", "D", "E", "F"), ("L", "M", "N", "O"), ("U", "V", "W", "X")]
        rep_tots, pcts = [], []
        for blk in blocks:
            c_cat, c_sum = blk[cat - 1], blk[3]
            img_pcts, cat_sum, tot_sum = [], 0.0, 0.0
            for r in range(r0, r1 + 1):
                if r in rows and c_sum in rows[r]:
                    s = float(rows[r][c_sum])
                    v = float(rows[r][c_cat])
                    cat_sum += v
                    tot_sum += s
                    img_pcts.append(100.0 * v / s)
            rep_tots.append(tot_sum)
            pcts.append(statistics.mean(img_pcts) if method == "image_mean" else 100.0 * cat_sum / tot_sum)

    elif ds.startswith("fig6d_"):
        rows = _read_sheet(src / "41567_2022_1924_MOESM17_ESM.xlsx", "Fig 6d")
        r0 = 7 if ds == "fig6d_3peg" else 20
        r_cat = r0 + (cat - 1)
        r_tot = r0 + 4
        img_blocks = [("C", "D", "E"), ("I", "J", "K"), ("O", "P", "Q")]
        sum_cols = ["F", "L", "R"]
        rep_tots = [float(rows[r_tot][c]) for c in sum_cols]
        if method == "image_mean":
            pcts = [
                statistics.mean(100.0 * float(rows[r_cat][c]) / float(rows[r_tot][c]) for c in blk)
                for blk in img_blocks
            ]
        else:
            pcts = [100.0 * float(rows[r_cat][c]) / float(rows[r_tot][c]) for c in sum_cols]

    elif ds.startswith("edfig4c_"):
        rows = _read_sheet(src / "41567_2022_1924_MOESM20_ESM.xlsx", "Fig S4c")
        cfg = {
            "edfig4c_1to20": ([6, 7, 8], 10),
            "edfig4c_1to10": ([15, 16, 17], 19),
            "edfig4c_1to7p5": ([24, 25, 26], 28),
        }
        img_rows, sum_row = cfg[ds]
        blocks = [("B", "C", "D", "E"), ("H", "I", "J", "K"), ("N", "O", "P", "Q")]
        rep_tots, pcts = [], []
        for blk in blocks:
            c_cat, c_sum = blk[cat - 1], blk[3]
            rep_tots.append(float(rows[sum_row][c_sum]))
            if method == "image_mean":
                pcts.append(statistics.mean(100.0 * float(rows[r][c_cat]) / float(rows[r][c_sum]) for r in img_rows))
            else:
                pcts.append(100.0 * float(rows[sum_row][c_cat]) / float(rows[sum_row][c_sum]))
    else:
        raise ValueError(f"unknown dataset: {ds!r}")

    return {
        "abundance_mean": statistics.mean(pcts),
        "abundance_sd": statistics.stdev(pcts),
        "total_droplets": sum(rep_tots),
        "min_replicate_droplets": min(rep_tots),
        "max_replicate_droplets": max(rep_tots),
    }
