"""Reads FigureS10C, FigureS10D, and FigureS10E (=published Supplementary Fig. 6C-E)."""
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


def _invert_cdf(xs: list[float], ys: list[float], q: float) -> float:
    for i in range(len(ys) - 1):
        if ys[i] <= q <= ys[i + 1] and ys[i + 1] > ys[i]:
            return xs[i] + (q - ys[i]) / (ys[i + 1] - ys[i]) * (xs[i + 1] - xs[i])
    return xs[-1]


@functools.lru_cache(maxsize=2)
def _load_nuc(src_str: str):
    sh = _read_sheets(Path(src_str) / _XLSX, ("FigureS10C", "FigureS10D", "FigureS10E"))
    max_int = 0.0

    # 1. FigureS10C filament length CDFs
    rs10c = sh["FigureS10C"]
    cdfs: dict[float, tuple[list[float], list[float]]] = {}
    cur = None
    for rn in range(2, 802):
        a = str(rs10c[rn]["A"]).strip()
        if a.startswith("[Nucleator]"):
            cur = float(a.split("=")[1].replace("nM", "").strip())
            cdfs[cur] = ([], [])
        else:
            cdfs[cur][0].append(_f(rs10c[rn]["A"]))
            cdfs[cur][1].append(_f(rs10c[rn]["B"]))

    # 2. FigureS10D VASP valency (93 raw columns = 3 reps x 31 frames, contradicting Supp. Fig. 6 caption '5 replicates x 30')
    rs10d = sh["FigureS10D"]
    vals_d: dict[float, dict[int, float]] = defaultdict(dict)
    cur = None
    n_raw_d = 0
    for rn in sorted(rs10d):
        a = str(rs10d[rn].get("A", "")).strip()
        if a.startswith("[Nucleator]"):
            cur = float(a.split("=")[1].replace("nM", "").strip())
        elif cur is not None and a in ("0", "1", "2", "3", "4"):
            m_st, s_st = _f(rs10d[rn]["B"]), _f(rs10d[rn]["C"])
            cols = sorted(
                [c for c in rs10d[rn] if c not in ("A", "B", "C") and str(rs10d[rn][c]).strip() != ""],
                key=_col_idx,
            )
            raw = [_f(rs10d[rn][c]) for c in cols]
            n_raw_d = len(raw)
            m_c, s_c = statistics.mean(raw), statistics.pstdev(raw)
            max_int = max(max_int, abs(m_c - m_st), abs(s_c - s_st))
            vals_d[cur][int(a)] = m_c

    # 3. FigureS10E surface area fraction (93 raw columns)
    rs10e = sh["FigureS10E"]
    vals_e: dict[float, tuple[float, float, int]] = {}
    for rn in range(2, 10):
        conc = _f(rs10e[rn]["A"])
        m_st, s_st = _f(rs10e[rn]["B"]), _f(rs10e[rn]["C"])
        cols = sorted(
            [c for c in rs10e[rn] if c not in ("A", "B", "C") and str(rs10e[rn][c]).strip() != ""],
            key=_col_idx,
        )
        raw = [_f(rs10e[rn][c]) for c in cols]
        m_c, s_c = statistics.mean(raw), statistics.pstdev(raw)
        max_int = max(max_int, abs(m_c - m_st), abs(s_c - s_st))
        vals_e[conc] = (m_c, s_c, len(raw))

    return cdfs, vals_d, vals_e, n_raw_d, max_int


def run(inputs, settings):
    cdfs, vals_d, vals_e, n_raw_d, max_int = _load_nuc(os.environ["BOXGRAPH_SOURCE"])
    cq = float(inputs["nuc_conc"])
    key = min(vals_e, key=lambda k: abs(k - cq))
    xs, ys = cdfs[key]
    first_1 = next((x for x, y in zip(xs, ys) if y >= 0.9999), xs[-1])
    af_m, af_s, n_raw_e = vals_e[key]
    v = vals_d[key]
    return {
        "length_cdf96": _invert_cdf(xs, ys, 0.96),
        "length_median": _invert_cdf(xs, ys, 0.50),
        "length_first_cdf1": first_1,
        "frac_0": v[0],
        "frac_1": v[1],
        "frac_2": v[2],
        "frac_3": v[3],
        "frac_4": v[4],
        "surface_area_mean": af_m,
        "surface_area_sd_pop": af_s,
        "n_raw_snapshots": float(n_raw_e),
        "internal_raw_vs_stored_max_diff": max_int,
    }
