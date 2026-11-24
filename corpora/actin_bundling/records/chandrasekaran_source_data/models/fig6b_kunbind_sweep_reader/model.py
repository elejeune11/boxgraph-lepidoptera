"""Reads sheet 'Figure6B' directly from 41467_2024_46726_MOESM9_ESM.xlsx in $BOXGRAPH_SOURCE."""
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
def _load_6b(src_str: str):
    rows = _read_sheet(Path(src_str) / _XLSX, "Figure6B")
    row5_raw = _f(rows[5]["E"])
    entries = []
    for rn in range(2, 16):
        ku_raw = _f(rows[rn]["E"])
        # Cell E5 misprints 1.15625 as 1.5625 (same typo printed on Fig. 6B x-axis); map by row position for lookup
        ku_pos = 1.15625 if rn == 5 and abs(ku_raw - 1.5625) < 1e-6 else ku_raw
        entries.append((ku_pos, 100.0 * _f(rows[rn]["D"]), 100.0 * _f(rows[rn]["C"]), 100.0 * _f(rows[rn]["B"])))
    return entries, row5_raw


def run(inputs, settings):
    entries, row5_raw = _load_6b(os.environ["BOXGRAPH_SOURCE"])
    ku_q = float(inputs["k_unbind"])
    _, rng, ssb, sh = min(entries, key=lambda e: abs(e[0] - ku_q))
    return {
        "ring_pct": rng,
        "shell_sb_pct": ssb,
        "shell_pct": sh,
        "n_rows": float(len(entries)),
        "row5_stored_kunbind": row5_raw,
    }
