"""Reads FigureS4A/B, FigureS8A/B (=Supp. Fig. 9A/B), FigureS11Bi/Bii, and FigureS12A/B."""
from __future__ import annotations

import functools
import os
import zipfile
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


def _varimax_communality_pct(rows: dict[int, dict[str, object]]) -> float:
    sq_sum = 0.0
    for rn in sorted(rows):
        if str(rows[rn].get("A", "")).strip().startswith("PC"):
            for col in ("B", "C", "D", "E", "F"):
                v = _f(rows[rn][col])
                sq_sum += v * v
    return 100.0 * sq_sum / 5.0


@functools.lru_cache(maxsize=2)
def _load_pca(src_str: str):
    names = ("FigureS4A", "FigureS4B", "FigureS8A", "FigureS8B", "FigureS11Bi", "FigureS11Bii", "FigureS12A", "FigureS12B")
    sh = _read_sheets(Path(src_str) / _XLSX, names)

    def _cum(sheet_name: str, cols: tuple[str, ...]) -> list[float]:
        r2 = sh[sheet_name][2]
        vals = [_f(r2[c]) for c in cols]
        return [sum(vals[: i + 1]) for i in range(len(vals))]

    cum_map = {
        "s4": _cum("FigureS4A", ("A", "B", "C", "D", "E")),
        "s9": _cum("FigureS8A", ("A", "B", "C", "D", "E")),
        "s12": _cum("FigureS12A", ("B", "C", "D", "E", "F")),
    }
    comm_map = {
        "s4": _varimax_communality_pct(sh["FigureS4B"]),
        "s9": _varimax_communality_pct(sh["FigureS8B"]),
        "s11": _varimax_communality_pct(sh["FigureS11Bii"]),
        "s12": _varimax_communality_pct(sh["FigureS12B"]),
    }

    # Verify FigureS11Bi == FigureS11Bii cell-for-cell (FigureS11Bi was overwritten with Varimax matrix!)
    s11_diff = 0.0
    for rn in (3, 4, 5):
        for col in ("B", "C", "D", "E", "F"):
            s11_diff = max(s11_diff, abs(_f(sh["FigureS11Bi"][rn][col]) - _f(sh["FigureS11Bii"][rn][col])))
    return cum_map, comm_map, s11_diff


def run(inputs, settings):
    cum_map, comm_map, s11_diff = _load_pca(os.environ["BOXGRAPH_SOURCE"])
    study = str(settings.get("study", "s4"))
    idx = int(round(float(inputs["pc_index"])))
    if study not in cum_map:
        raise ValueError(f"Unsupported study {study!r}; use 's4', 's9', or 's12'")
    cum_r = cum_map[study][idx - 1]
    return {
        "cum_var_ratio": cum_r,
        "cum_var_pct": 100.0 * cum_r,
        "varimax_mean_communality_pct": comm_map[study],
        "s11_varimax_mean_communality_pct": comm_map["s11"],
        "s11bi_vs_s11bii_max_diff": s11_diff,
    }
