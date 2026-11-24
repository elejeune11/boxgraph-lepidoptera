"""Interpolates live average FRAP recovery curves from MOESM12 ('Fig 1f') and MOESM17 ('Fig 6b')."""
from __future__ import annotations

import os
import statistics
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np

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
    ds = settings.get("dataset", "fig1f")
    if ds == "fig1f":
        rows = _read_sheet(src / "41567_2022_1924_MOESM12_ESM.xlsx", "Fig 1f")
        tcol, cols, r0 = "B", ["C", "D", "E", "F", "G", "H", "I", "J", "K"], 4
    elif ds == "fig6b_3peg":
        rows = _read_sheet(src / "41567_2022_1924_MOESM17_ESM.xlsx", "Fig 6b")
        tcol, cols, r0 = "A", ["B", "C", "D", "E", "F"], 5
    elif ds == "fig6b_10peg":
        rows = _read_sheet(src / "41567_2022_1924_MOESM17_ESM.xlsx", "Fig 6b")
        tcol, cols, r0 = "H", ["I", "J", "K", "L", "M"], 5
    else:
        raise ValueError(f"unknown dataset: {ds!r}")

    ts, ys = [], []
    for r in sorted(rows):
        if r >= r0 and tcol in rows[r]:
            vs = [float(rows[r][c]) * 100.0 for c in cols if c in rows[r]]
            if vs:
                ts.append(float(rows[r][tcol]))
                ys.append(statistics.mean(vs))
    t_query = float(inputs["t"])
    return {"recovery": float(np.interp(t_query, ts, ys))}
