"""Live software adapter reading PLOTCSVFILES/Fillength_VCA_Arp23_concsweep_steric_kunbind_0R_*.csv."""

import csv
import math
import os
from pathlib import Path
import scipy.constants

_ARP_CONCS = [0.0, 10.0, 25.0, 50.0, 75.0, 100.0, 125.0, 150.0]
_VCA_CONCS = [0.0, 1.0, 2.0, 5.0, 10.0]


def _find_idx(val, grid, name):
    for idx, g in enumerate(grid):
        if abs(val - g) < 1e-6:
            return idx
    raise ValueError(f"{name}={val} not in grid {grid}")


def _read_fillength(csv_dir, rval):
    fp = csv_dir / f"Fillength_VCA_Arp23_concsweep_steric_kunbind_0R_{rval}.csv"
    with fp.open(encoding="utf-8") as fh:
        return {row[0]: [float(x) for x in row[1:]] for row in csv.reader(fh) if row and row[0]}


def _infer_total_vca_fillength_reps(csv_dir):
    # Check denominator of Nfil_mean fractional parts across all 40 Fillength_VCA CSVs:
    # if 3 replicates were averaged, multiples of 1/3 (~0.3333, 0.6667) must appear;
    # if only 2 replicates were averaged, only {0.0, 0.5} appear.
    fracs = set()
    for rval in range(40):
        rows = _read_fillength(csv_dir, rval)
        for v in rows["Nfil_mean"]:
            fracs.add(round(v - math.floor(v), 4))
    if any(abs(f - 1.0 / 3.0) < 0.02 or abs(f - 2.0 / 3.0) < 0.02 for f in fracs):
        return 3.0
    if 0.5 in fracs:
        return 2.0
    return 1.0


def run(inputs, settings):
    a_idx = _find_idx(float(inputs["C_Arp"]), _ARP_CONCS, "C_Arp")
    v_idx = _find_idx(float(inputs["C_VCA"]), _VCA_CONCS, "C_VCA")
    R = float(inputs["R_drop"])
    l_mono = float(inputs["l_mono"])

    rval = a_idx * 5 + v_idx
    rval_base = a_idx * 5  # [VCA] = 0 uM at same [Arp2/3]
    csv_dir = Path(os.environ["BOXGRAPH_SOURCE"]) / "PLOTCSVFILES"
    rows = _read_fillength(csv_dir, rval)
    rows_base = _read_fillength(csv_dir, rval_base)

    nfil_0 = rows["Nfil_mean"][0]
    nfil_600 = rows["Nfil_mean"][-1]
    lfil_600 = rows["Lfil_mean"][-1]
    lsum_600 = rows["Lsum_mean"][-1]
    lsum_base_600 = rows_base["Lsum_mean"][-1]

    vol_L = (4.0 / 3.0) * math.pi * (R ** 3) * 1e-15
    uM_per_um = (1.0 / l_mono) / (scipy.constants.Avogadro * vol_L) * 1e6

    inferred_reps = _infer_total_vca_fillength_reps(csv_dir)

    return {
        "N_fil_600": float(nfil_600),
        "N_fil_fold": float(nfil_600 / nfil_0),
        "L_fil_mean_600": float(lfil_600),
        "L_sum_600": float(lsum_600),
        "F_actin_conc_600": float(lsum_600 * uM_per_um),
        "delta_F_actin_600": float((lsum_base_600 - lsum_600) * uM_per_um),
        "inferred_csv_replicates": float(inferred_reps),
    }
