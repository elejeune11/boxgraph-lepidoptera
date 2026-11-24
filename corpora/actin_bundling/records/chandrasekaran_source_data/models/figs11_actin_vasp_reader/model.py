"""Reads FigureS11C and cross-checks against FigureS11 (=published Supplementary Fig. 11D)."""
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
def _load_s11(src_str: str):
    sh = _read_sheets(Path(src_str) / _XLSX, ("FigureS11C", "FigureS11"))
    r11c = sh["FigureS11C"]
    counts: dict[tuple[float, int], Counter] = defaultdict(Counter)
    reps: dict[tuple[float, int], set[int]] = defaultdict(set)
    for rn in range(7, 497):
        cl = int(_f(r11c[rn]["C"]))
        tag = ast.literal_eval(str(r11c[rn]["D"]))
        if tag[0] == "Actin:VASP":
            _, vasp, nfil, rep, _ = tag
            counts[(float(vasp), int(nfil))][cl] += 1
            reps[(float(vasp), int(nfil))].add(int(rep))

    r11 = sh["FigureS11"]
    cur_vasp = None
    max_diff = 0.0
    for rn in range(1, 20):
        a = str(r11[rn].get("A", "")).strip()
        if a.startswith("[VASP]="):
            cur_vasp = float(a.split("=")[1])
        elif "B" in r11[rn] and str(r11[rn]["B"]).strip() not in ("", "Column2", "Nfil"):
            nfil = int(round(_f(r11[rn]["B"])))
            st = [_f(r11[rn]["C"]), _f(r11[rn]["D"]), _f(r11[rn]["E"])]
            c = counts[(cur_vasp, nfil)]
            tot = sum(c.values())
            comp = [c[k] / tot for k in (0, 1, 2)]
            max_diff = max(max_diff, max(abs(x - y) for x, y in zip(st, comp)))
    return counts, reps, max_diff


def run(inputs, settings):
    counts, reps, max_diff = _load_s11(os.environ["BOXGRAPH_SOURCE"])
    vasp = float(settings.get("vasp_uM", 2.0))
    nfil = int(round(float(inputs["n_fil"])))
    key = min(counts, key=lambda k: abs(k[0] - vasp) * 1000.0 + abs(k[1] - nfil))
    c = counts[key]
    tot = sum(c.values())
    return {
        "ring_shell_sb_prob": c[0] / tot,
        "ring_multiring_prob": c[1] / tot,
        "locally_bundled_prob": c[2] / tot,
        "n_snapshots": float(tot),
        "n_replicates": float(len(reps[key])),
        "s11c_vs_s11_max_diff": max_diff,
    }
