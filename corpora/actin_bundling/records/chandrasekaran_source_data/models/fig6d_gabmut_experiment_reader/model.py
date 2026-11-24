"""Reads sheet 'Figure6D' (VASP-GABmut experimental replicates 1-3; WT rows 1-19 are blank)."""
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
def _load_6d(src_str: str):
    rows = _read_sheet(Path(src_str) / _XLSX, "Figure6D")
    # Three Actin:VASP ratios: 1:20 (ratio=20), 1:10 (ratio=10), 1:7.5 (ratio=7.5)
    # Replicate count rows: rep1=(23,24,25), rep2=(32,33,34), rep3=(41,42,43)
    # Stored summary rows: AVERAGE=(22,23,24), SD=(30,31,32), SEM=(36,37,38) in cols L,M,N
    specs = [
        (20.0, [23, 32, 41], 22, 30, 36),
        (10.0, [24, 33, 42], 23, 31, 37),
        (7.5,  [25, 34, 43], 24, 32, 38),
    ]
    data = {}
    max_msd_diff = 0.0
    max_sem_div3_diff = 0.0
    total_all = 0.0
    for ratio, rep_rows, avg_r, sd_r, sem_r in specs:
        sh_pcts, rg_pcts, dk_pcts, ns = [], [], [], []
        for rr in rep_rows:
            sh, rg, dk = _f(rows[rr]["B"]), _f(rows[rr]["C"]), _f(rows[rr]["D"])
            tot = sh + rg + dk
            ns.append(tot)
            sh_pcts.append(100.0 * sh / tot)
            rg_pcts.append(100.0 * rg / tot)
            dk_pcts.append(100.0 * dk / tot)
            max_msd_diff = max(
                max_msd_diff,
                abs(tot - _f(rows[rr]["E"])),
                abs(sh_pcts[-1] - _f(rows[rr]["F"])),
                abs(rg_pcts[-1] - _f(rows[rr]["G"])),
                abs(dk_pcts[-1] - _f(rows[rr]["H"])),
            )
        total_all += sum(ns)
        m_sh, m_rg, m_dk = statistics.mean(sh_pcts), statistics.mean(rg_pcts), statistics.mean(dk_pcts)
        s_sh, s_rg, s_dk = statistics.stdev(sh_pcts), statistics.stdev(rg_pcts), statistics.stdev(dk_pcts)
        st_sem = (_f(rows[sem_r]["L"]), _f(rows[sem_r]["M"]), _f(rows[sem_r]["N"]))
        max_msd_diff = max(
            max_msd_diff,
            abs(m_sh - _f(rows[avg_r]["L"])),
            abs(m_rg - _f(rows[avg_r]["M"])),
            abs(m_dk - _f(rows[avg_r]["N"])),
            abs(s_sh - _f(rows[sd_r]["L"])),
            abs(s_rg - _f(rows[sd_r]["M"])),
            abs(s_dk - _f(rows[sd_r]["N"])),
        )
        max_sem_div3_diff = max(
            max_sem_div3_diff,
            abs(s_sh / 3.0 - st_sem[0]),
            abs(s_rg / 3.0 - st_sem[1]),
            abs(s_dk / 3.0 - st_sem[2]),
        )
        sq3 = math.sqrt(len(rep_rows))
        data[ratio] = {
            "shell_mean": m_sh,
            "ring_mean": m_rg,
            "disk_mean": m_dk,
            "shell_sd": s_sh,
            "ring_sd": s_rg,
            "disk_sd": s_dk,
            "shell_sem_true": s_sh / sq3,
            "ring_sem_true": s_rg / sq3,
            "disk_sem_true": s_dk / sq3,
            "shell_sem_stored": st_sem[0],
            "ring_sem_stored": st_sem[1],
            "disk_sem_stored": st_sem[2],
            "n_replicates": float(len(rep_rows)),
            "n_droplets_condition": sum(ns),
        }
    return data, total_all, max_msd_diff, max_sem_div3_diff


def run(inputs, settings):
    data, total_all, max_msd_diff, max_sem_div3_diff = _load_6d(os.environ["BOXGRAPH_SOURCE"])
    r_q = float(inputs["vasp_to_actin_ratio"])
    key = min(data, key=lambda k: abs(k - r_q))
    out = dict(data[key])
    out["n_droplets_all_conditions"] = total_all
    out["internal_mean_sd_max_diff"] = max_msd_diff
    out["stored_sem_vs_sd_over_3_max_diff"] = max_sem_div3_diff
    return out
