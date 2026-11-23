"""Reads Source-data_SupplementaryFig3a.xlsx directly from $BOXGRAPH_SOURCE."""
from __future__ import annotations

import os
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

_M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


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


def run(inputs, settings):
    src = Path(os.environ["BOXGRAPH_SOURCE"])
    genotype = str(settings.get("genotype", "attP40"))
    rows = _read_sheet(src / "Source-data_SupplementaryFig3a.xlsx", "Supplementary Fig 3a")

    headers = {str(v).strip(): k for k, v in rows[7].items() if k != "B"}
    if genotype not in headers:
        raise KeyError(f"genotype {genotype!r} not in SupplementaryFig3a ({sorted(headers)})")
    col = headers[genotype]

    pairs = sorted((float(rows[r]["B"]), float(rows[r][col])) for r in range(8, 98))
    tot = sum(y for _, y in pairs)
    cum = 0.0
    med_asc = None
    med_prev = None
    med_interp = None
    for a, y in pairs:
        prev_cum = cum
        cum += y
        if cum >= 0.5 * tot and med_asc is None:
            med_asc = a
            med_prev = a - 1.0
            # Bin [a - 0.5, a + 0.5]
            frac = (0.5 * tot - prev_cum) / y if y > 0 else 0.5
            med_interp = (a - 0.5) + frac
    return {
        "median_orientation_asc": float(med_asc),
        "median_orientation_prev": float(med_prev),
        "median_orientation_interp": float(med_interp),
        "pdf_sum": tot,
        "n_genotypes_present": float(len(headers)),
    }
