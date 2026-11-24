"""Reads Figure5C, Figure5D, and FigureS6 directly from 41467_2024_46726_MOESM9_ESM.xlsx."""
from __future__ import annotations

import ast
import functools
import os
import zipfile
from collections import Counter, defaultdict
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
def _load_fig5(src_str: str):
    sh = _read_sheets(Path(src_str) / _XLSX, ("Figure5C", "Figure5D", "FigureS6"))
    r5c = sh["Figure5C"]
    counts: dict[tuple[float, float], Counter] = defaultdict(Counter)
    reps_by_cell: dict[tuple[float, float], set[int]] = defaultdict(set)
    for rn in range(2, len(r5c) + 1):
        cl = int(_f(r5c[rn]["E"]))
        tag = ast.literal_eval(str(r5c[rn]["F"]))
        if tag[0] == "T2piR":
            _, vasp, t2p, rep, _ = tag
            counts[(float(vasp), float(t2p))][cl] += 1
            reps_by_cell[(float(vasp), float(t2p))].add(int(rep))

    # In Figure5C, Cluster 1 = Ring and Cluster 0 = Shell-SB (matching 30/30 Figure5D entries)
    r5d = sh["Figure5D"]
    t2p_cols = [(10.0, "B"), (50.0, "C"), (150.0, "D"), (200.0, "E"), (450.0, "F"), (600.0, "G")]
    stored_5d: dict[tuple[float, float], float] = {}
    max_diff = 0.0
    for rn in range(3, 8):
        vasp = _f(r5d[rn]["A"])
        for t2p, col in t2p_cols:
            p5d = 100.0 * _f(r5d[rn][col])
            stored_5d[(vasp, t2p)] = p5d
            c = counts[(vasp, t2p)]
            p5c = 100.0 * c[1] / sum(c.values())
            max_diff = max(max_diff, abs(p5c - p5d))

    rs6 = sh["FigureS6"]
    s6_nonzero = sum(1 for rn in range(2, len(rs6) + 1) if int(_f(rs6[rn]["D"])) != 0)
    return counts, reps_by_cell, stored_5d, max_diff, s6_nonzero


def run(inputs, settings):
    counts, reps_by_cell, stored_5d, max_diff, s6_nonzero = _load_fig5(os.environ["BOXGRAPH_SOURCE"])
    vasp = float(settings.get("vasp_uM", 0.8))
    t2p = float(inputs["t_2pir"])
    key = min(stored_5d, key=lambda k: abs(k[0] - vasp) * 1000.0 + abs(k[1] - t2p))
    c = counts[key]
    tot = sum(c.values())
    return {
        "ring_prob_from_5c": 100.0 * c[1] / tot,
        "ring_prob_stored_5d": stored_5d[key],
        "n_snapshots_per_cell": float(tot),
        "n_replicates": float(len(reps_by_cell[key])),
        "max_5c_vs_5d_diff": max_diff,
        "figs6_nonzero_cluster_count": float(s6_nonzero),
    }
