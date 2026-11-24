"""Reads FigureS5 and FigureS9 (=Supp. Fig. 10), decoding inverted span headers (spanc>=spanb>=spana)
and scrambled workbook cluster legends."""
from __future__ import annotations

import functools
import os
import statistics
import zipfile
from collections import defaultdict
from pathlib import Path
from xml.etree import ElementTree as ET

_M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
_XLSX = "41467_2024_46726_MOESM9_ESM.xlsx"


def _f(x: object) -> float:
    return float(str(x).strip())


def _read_sheets(path: Path, names: tuple[str, ...]) -> dict[str, dict[int, dict[str, object]]]:
    out: dict[str, dict[int, dict[str, object]]] = {}
    with zipfile.ZipFile(path) as z:
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        relmap = {r.get("Id"): r.get("Target") for r in rels}
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            for si in ET.fromstring(z.read("xl/sharedStrings.xml")).iter(_M + "si"):
                shared.append("".join(t.text or "" for t in si.iter(_M + "t")))
        sheet_map = {s.get("name"): s.get(_R + "id") for s in wb.iter(_M + "sheet")}
        for name in names:
            target = relmap[sheet_map[name]].lstrip("/")
            if not target.startswith("xl/"):
                target = "xl/" + target
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
            out[name] = rows
    return out


@functools.lru_cache(maxsize=2)
def _load_spans(src_str: str):
    sh = _read_sheets(Path(src_str) / _XLSX, ("FigureS5", "FigureS9"))
    # FigureS5: col A='spanc' (largest -> plotted a), col B='spanb' (-> b), col C='spana' (smallest -> c)
    # Decoded clusters: 0=Shell, 3=Shell-WB, 1=Shell-SB, 2=Ring (Supp. Fig. 5 x-axis categories 1..4)
    s5_raw = defaultdict(lambda: ([], [], []))
    for rn in range(5, 545):
        cl = int(_f(sh["FigureS5"][rn]["D"]))
        s5_raw[cl][0].append(_f(sh["FigureS5"][rn]["A"]))
        s5_raw[cl][1].append(_f(sh["FigureS5"][rn]["B"]))
        s5_raw[cl][2].append(_f(sh["FigureS5"][rn]["C"]))
    s5_map = {1: 0, 2: 3, 3: 1, 4: 2}
    s5 = {
        cat: (
            statistics.median(s5_raw[cl][0]),
            statistics.median(s5_raw[cl][1]),
            statistics.median(s5_raw[cl][2]),
            len(s5_raw[cl][0]),
        )
        for cat, cl in s5_map.items()
    }

    # FigureS9 (=Supp. Fig. 10): decoded clusters: 1=Ring (flat c=0.144 um), 0=Shell-SB (c=0.736 um)
    # (Row 2 legend '0-Ring, 1-Shell-SB' in workbook is inverted!)
    s9_raw = defaultdict(lambda: ([], [], []))
    for rn in range(4, 1039):
        cl = int(_f(sh["FigureS9"][rn]["D"]))
        s9_raw[cl][0].append(_f(sh["FigureS9"][rn]["A"]))
        s9_raw[cl][1].append(_f(sh["FigureS9"][rn]["B"]))
        s9_raw[cl][2].append(_f(sh["FigureS9"][rn]["C"]))
    s10_map = {1: 1, 2: 0}
    s10 = {
        cat: (
            statistics.median(s9_raw[cl][0]),
            statistics.median(s9_raw[cl][1]),
            statistics.median(s9_raw[cl][2]),
            len(s9_raw[cl][0]),
        )
        for cat, cl in s10_map.items()
    }
    return {"s5": s5, "s10": s10}


def run(inputs, settings):
    data = _load_spans(os.environ["BOXGRAPH_SOURCE"])
    study = str(settings.get("study", "s5"))
    cat = int(round(float(inputs["category"])))
    a_m, b_m, c_m, n = data[study][cat]
    return {
        "span_a_median": a_m,
        "span_b_median": b_m,
        "span_c_median": c_m,
        "n_snapshots": float(n),
    }
