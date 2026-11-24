"""Reads VASP and actin partitioning tables from MOESM18 ('Fig S2d') and MOESM19 ('Fig S3e')."""
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
    ds = settings.get("dataset", "edfig2d")
    idx = int(round(float(inputs["construct"])))

    if ds == "edfig2d":
        rows = _read_sheet(src / "41567_2022_1924_MOESM18_ESM.xlsx", "Fig S2d")
        cols_map = {
            1: ["B", "C", "D"],
            2: ["F", "G", "H"],
            3: ["J", "K", "L"],
            4: ["N", "O", "P"],
        }
        cols = cols_map[idx]
    elif ds == "edfig3e":
        rows = _read_sheet(src / "41567_2022_1924_MOESM19_ESM.xlsx", "Fig S3e")
        cols_map = {
            1: ["C", "D", "E"],
            2: ["G", "H", "I"],
        }
        cols = cols_map[idx]
    else:
        raise ValueError(f"unknown dataset: {ds!r}")

    trial_means, all_vals = [], []
    for c in cols:
        vs = [float(rows[r][c]) for r in sorted(rows) if r >= 3 and c in rows[r]]
        trial_means.append(statistics.mean(vs))
        all_vals.extend(vs)

    g_mean = statistics.mean(trial_means)
    g_sd = statistics.stdev(trial_means)
    g_sem = g_sd / math.sqrt(len(trial_means))
    return {
        "partition_mean": g_mean,
        "partition_sd": g_sd,
        "partition_sem": g_sem,
        "pooled_mean": statistics.mean(all_vals),
        "n_droplets_total": float(len(all_vals)),
    }
