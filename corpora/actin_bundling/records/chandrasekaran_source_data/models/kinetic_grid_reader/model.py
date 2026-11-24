"""Reads FigureS2C, Figure3D, Figure3E, and cross-checks FigureS2C against Figure3B (t=570..600 s)."""
from __future__ import annotations

import ast
import functools
import math
import os
import statistics
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from xml.etree import ElementTree as ET

_M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
_XLSX = "41467_2024_46726_MOESM9_ESM.xlsx"


def _col_idx(col: str) -> int:
    n = 0
    for ch in col:
        n = n * 26 + (ord(ch) - 64)
    return n


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
def _load_grid(src_str: str):
    sh = _read_sheets(Path(src_str) / _XLSX, ("Figure3B", "FigureS2C", "Figure3D", "Figure3E"))

    # 1. Figure3B last 31 frames (t=570..600 s inclusive) x 3 replicates = 93 values
    r3b = sh["Figure3B"]
    f3b_tail: dict[tuple[float, float], list[float]] = defaultdict(list)
    cur = None
    for rn in sorted(r3b):
        a = str(r3b[rn].get("A", "")).strip()
        if a.startswith("kbind="):
            cur = (float(a.split("=")[1]), float(str(r3b[rn]["B"]).strip().split("=")[1]))
        elif cur is not None and a not in ("Time", "Column1", ""):
            if _f(r3b[rn]["A"]) >= 570.0:
                f3b_tail[cur].extend([_f(r3b[rn]["B"]), _f(r3b[rn]["C"]), _f(r3b[rn]["D"])])

    # 2. FigureS2C raw 93 columns -> recompute mean, pop SD, and compare to Figure3B tail
    rs2c = sh["FigureS2C"]
    s2c: dict[tuple[float, float], tuple[float, float, int]] = {}
    max_s2c_3b = 0.0
    max_int = 0.0
    for rn in range(2, 38):
        pair = (_f(rs2c[rn]["A"]), _f(rs2c[rn]["B"]))
        m_st, s_st = _f(rs2c[rn]["C"]), _f(rs2c[rn]["D"])
        cols = sorted(
            [c for c in rs2c[rn] if c not in ("A", "B", "C", "D") and str(rs2c[rn][c]).strip() != ""],
            key=_col_idx,
        )
        raw = [_f(rs2c[rn][c]) for c in cols]
        m_c, s_c = statistics.mean(raw), statistics.pstdev(raw)
        max_int = max(max_int, abs(m_c - m_st), abs(s_c - s_st))
        tail = sorted(f3b_tail[pair])
        if len(tail) == len(raw):
            max_s2c_3b = max(max_s2c_3b, max(abs(x - y) for x, y in zip(tail, sorted(raw))))
        s2c[pair] = (m_c, s_c, len(raw))

    # 3. Figure3D raw 93 columns -> recompute valency 0..4 means & pop SDs
    r3d = sh["Figure3D"]
    f3d: dict[tuple[float, float], dict[int, tuple[float, float]]] = defaultdict(dict)
    cur = None
    for rn in sorted(r3d):
        a = str(r3d[rn].get("A", "")).strip()
        if a.startswith("kbind="):
            cur = (float(a.split("=")[1]), float(str(r3d[rn]["B"]).strip().split("=")[1]))
        elif cur is not None and a in ("0", "1", "2", "3", "4"):
            valency = int(a)
            m_st, s_st = _f(r3d[rn]["B"]), _f(r3d[rn]["C"])
            cols = sorted(
                [c for c in r3d[rn] if c not in ("A", "B", "C") and str(r3d[rn][c]).strip() != ""],
                key=_col_idx,
            )
            raw = [_f(r3d[rn][c]) for c in cols]
            m_c, s_c = statistics.mean(raw), statistics.pstdev(raw)
            max_int = max(max_int, abs(m_c - m_st), abs(s_c - s_st))
            f3d[cur][valency] = (m_c, s_c)

    # 4. Figure3E cluster IDs (actual decoded cluster ids: 0=Shell, 1=Shell-SB, 2=Ring, 3=Shell-WB;
    #    note Row 2 legend in the workbook is scrambled!)
    r3e = sh["Figure3E"]
    f3e: dict[tuple[float, float], Counter] = defaultdict(Counter)
    for rn in range(6, 546):
        cl = int(_f(r3e[rn]["D"]))
        kb, ku, _, _ = ast.literal_eval(str(r3e[rn]["E"]))
        f3e[(float(kb), float(ku))][cl] += 1

    return s2c, f3d, f3e, max_s2c_3b, max_int


def run(inputs, settings):
    s2c, f3d, f3e, max_s2c_3b, max_int = _load_grid(os.environ["BOXGRAPH_SOURCE"])
    kb = float(settings.get("k_bind", 1.0))
    ku = float(inputs["k_unbind"])
    pair = min(
        s2c,
        key=lambda k: abs(math.log10(k[0]) - math.log10(kb)) + abs(math.log10(k[1]) - math.log10(ku)),
    )
    af_m, af_s, n_raw = s2c[pair]
    v = f3d[pair]
    cnt = f3e[pair]
    tot = sum(cnt.values())
    return {
        "surface_area_mean": af_m,
        "surface_area_sd_pop": af_s,
        "frac_0": v[0][0],
        "frac_1": v[1][0],
        "frac_2": v[2][0],
        "frac_3": v[3][0],
        "frac_4": v[4][0],
        "valency_sd_0": v[0][1],
        "valency_sd_4": v[4][1],
        "shell_pct": 100.0 * cnt[0] / tot,
        "shell_wb_pct": 100.0 * cnt[3] / tot,
        "shell_sb_pct": 100.0 * cnt[1] / tot,
        "ring_pct": 100.0 * cnt[2] / tot,
        "n_raw_snapshots": float(n_raw),
        "n_cluster_snapshots": float(tot),
        "s2c_vs_fig3b_last31_max_diff": max_s2c_3b,
        "internal_raw_vs_stored_max_diff": max_int,
    }
