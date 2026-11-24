"""Reads sheet 'Figure3B' directly from 41467_2024_46726_MOESM9_ESM.xlsx in $BOXGRAPH_SOURCE."""
from __future__ import annotations

import functools
import os
import statistics
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

_M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
_XLSX = "41467_2024_46726_MOESM9_ESM.xlsx"


def _f(x: object) -> float:
    return float(str(x).strip())


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


@functools.lru_cache(maxsize=2)
def _load_blocks(src_str: str):
    rows = _read_sheet(Path(src_str) / _XLSX, "Figure3B")
    blocks: dict[tuple[float, float], list[tuple[float, float, float, float]]] = {}
    cur = None
    max_m_err = 0.0
    max_s_err = 0.0
    for rn in sorted(rows):
        a = str(rows[rn].get("A", "")).strip()
        if a.startswith("kbind="):
            kb = float(a.split("=")[1])
            ku = float(str(rows[rn]["B"]).strip().split("=")[1])
            cur = (kb, ku)
            blocks[cur] = []
        elif cur is not None and a not in ("Time", "Column1", ""):
            t = _f(rows[rn]["A"])
            reps = [_f(rows[rn]["B"]), _f(rows[rn]["C"]), _f(rows[rn]["D"])]
            m_st = _f(rows[rn]["E"])
            s_st = _f(rows[rn]["F"])
            m_c = statistics.mean(reps)
            s_c = statistics.pstdev(reps)
            max_m_err = max(max_m_err, abs(m_c - m_st))
            max_s_err = max(max_s_err, abs(s_c - s_st))
            blocks[cur].append((t, m_c, s_c, m_st))
    return blocks, max_m_err, max_s_err


def _interp(xs: list[float], ys: list[float], xq: float) -> float:
    if xq <= xs[0]:
        return ys[0]
    if xq >= xs[-1]:
        return ys[-1]
    for i in range(len(xs) - 1):
        if xs[i] <= xq <= xs[i + 1]:
            if xs[i + 1] == xs[i]:
                return ys[i]
            w = (xq - xs[i]) / (xs[i + 1] - xs[i])
            return ys[i] + w * (ys[i + 1] - ys[i])
    return ys[-1]


def run(inputs, settings):
    blocks, max_m_err, max_s_err = _load_blocks(os.environ["BOXGRAPH_SOURCE"])
    kb = float(settings.get("k_bind", 10.0))
    ku = float(settings.get("k_unbind", 1.0))
    key = min(blocks, key=lambda k: abs(k[0] - kb) + abs(k[1] - ku))
    series = blocks[key]
    ts = [p[0] for p in series]
    t_q = float(inputs["t"])
    return {
        "surface_area_mean": _interp(ts, [p[1] for p in series], t_q),
        "surface_area_sd_pop": _interp(ts, [p[2] for p in series], t_q),
        "surface_area_mean_stored": _interp(ts, [p[3] for p in series], t_q),
        "n_replicates": 3.0,
        "internal_mean_max_diff": max_m_err,
        "internal_sd_pop_max_diff": max_s_err,
    }
