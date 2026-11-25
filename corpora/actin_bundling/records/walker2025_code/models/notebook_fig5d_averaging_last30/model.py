"""Corrected last-30-snapshots variant of Figures.ipynb Cell 13 operating on $BOXGRAPH_SOURCE CSVs."""

import csv
import os
from pathlib import Path

import numpy as np


def run(inputs, settings):
    src_root = Path(os.environ["BOXGRAPH_SOURCE"])
    r_val = int(round(float(inputs["r_index"])))
    n_reps = int(round(float(inputs["n_reps"])))
    tv = float(inputs["n_monomers_total"])
    n_last = int(settings.get("n_last_snapshots", 30))

    # Row indices in Crosslink_inputfiles_R_*_r_*.csv (after header row):
    # 0=Time, 1=Bound_Dimer, 2=Bound_Monomer, 3=Free_Dimer, 4=Free_Monomer,
    # 5=Monomer_0, 6=Monomer_1, 7=Dimer_0, 8=Dimer_1, 9=Dimer_2
    buckets = {k: [] for k in ("frac_M0", "frac_M1", "frac_D0", "frac_D1", "frac_D2")}
    row_map = {"frac_M0": 5, "frac_M1": 6, "frac_D0": 7, "frac_D1": 8, "frac_D2": 9}

    for rep in range(n_reps):
        csv_path = (
            src_root
            / "Figure_5"
            / "outputfiles"
            / "CrosslinkOutput"
            / f"Crosslink_inputfiles_R_{r_val}_r_{rep}.csv"
        )
        with open(csv_path, newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            next(reader)
            rows = [np.array([float(x) for x in r[1:]], dtype=float) for r in reader]
        for key, ridx in row_map.items():
            buckets[key].append(rows[ridx][-n_last:] / tv)

    return {
        "n_snapshots_per_rep": float(n_last),
        "frac_M0": float(np.mean(np.concatenate(buckets["frac_M0"]))),
        "frac_M1": float(np.mean(np.concatenate(buckets["frac_M1"]))),
        "frac_D0": float(np.mean(np.concatenate(buckets["frac_D0"]))),
        "frac_D1": float(np.mean(np.concatenate(buckets["frac_D1"]))),
        "frac_D2": float(np.mean(np.concatenate(buckets["frac_D2"]))),
    }
