"""Reads FigureS4C, FigureS8C (=Supp. Fig. 9C), FigureS11Biii, and FigureS12C."""
from __future__ import annotations

import functools
import os
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

_M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
_XLSX = "41467_2024_46726_MOESM9_ESM.xlsx"
_SHEETS = {"s4": "FigureS4C", "s9": "FigureS8C", "s11": "FigureS11Biii", "s12": "FigureS12C"}


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
def _load_sil(src_str: str):
    sh = _read_sheets(Path(src_str) / _XLSX, tuple(_SHEETS.values()))
    out = {}
    for key, sname in _SHEETS.items():
        rows = sh[sname]
        table = {int(round(_f(rows[rn]["A"]))): _f(rows[rn]["B"]) for rn in range(2, 11)}
        opt_k = max(table, key=lambda k: table[k])
        out[key] = (table, opt_k, table[opt_k])
    return out


def run(inputs, settings):
    data = _load_sil(os.environ["BOXGRAPH_SOURCE"])
    study = str(settings.get("study", "s4"))
    table, opt_k, opt_s = data[study]
    k = int(round(float(inputs["k"])))
    return {
        "silhouette": table[k],
        "optimal_k": float(opt_k),
        "optimal_silhouette": opt_s,
    }
