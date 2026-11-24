"""Reads sheet 'Figure7B' (5 replicates per L_fil^max, confirming Methods p. 11 vs Fig. 7B caption's 'six replicates')."""
from __future__ import annotations

import functools
import math
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
def _load_7b(src_str: str):
    rows = _read_sheet(Path(src_str) / _XLSX, "Figure7B")
    entries = []
    max_diff = 0.0
    for rn in range(2, 12):
        L_um = _f(rows[rn]["A"]) * math.pi
        m_st, s_st = _f(rows[rn]["B"]), _f(rows[rn]["C"])
        raw = [_f(rows[rn][c]) for c in ("D", "E", "F", "G", "H")]
        m_c = statistics.mean(raw)
        s_p = statistics.pstdev(raw)
        s_s = statistics.stdev(raw)
        max_diff = max(max_diff, abs(m_c - m_st), abs(s_p - s_st))
        entries.append((L_um, m_c, s_p, s_s, raw))
    return entries, max_diff


def run(inputs, settings):
    entries, max_diff = _load_7b(os.environ["BOXGRAPH_SOURCE"])
    lq = float(inputs["l_fil_max"])
    _, m_c, s_p, s_s, raw = min(entries, key=lambda e: abs(e[0] - lq))
    return {
        "aspect_ratio_mean": m_c,
        "aspect_ratio_sd_pop": s_p,
        "aspect_ratio_sd_sample": s_s,
        "n_replicates": float(len(raw)),
        "rep_0": raw[0],
        "rep_1": raw[1],
        "rep_2": raw[2],
        "rep_3": raw[3],
        "rep_4": raw[4],
        "internal_mean_sd_pop_max_diff": max_diff,
    }
