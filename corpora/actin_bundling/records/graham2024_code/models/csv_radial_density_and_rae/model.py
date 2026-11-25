"""Live software adapter executing Cell 39 / Cell 43 radial density and R.A.E. logic directly on PLOTCSVFILES."""

import csv
import os
from pathlib import Path
import numpy as np

_ARP_CONCS = [0.0, 10.0, 25.0, 50.0, 75.0, 100.0, 125.0, 150.0]
_VCA_CONCS = [0.0, 1.0, 2.0, 5.0, 10.0]


def _find_idx(val, grid, name):
    for idx, g in enumerate(grid):
        if abs(val - g) < 1e-6:
            return idx
    raise ValueError(f"{name}={val} not in grid {grid}")


def run(inputs, settings):
    dataset = str(settings.get("dataset", "arp23"))
    n_tail = int(settings.get("n_tail_frames", 10))
    a_idx = _find_idx(float(inputs["C_Arp"]), _ARP_CONCS, "C_Arp")
    v_idx = _find_idx(float(inputs["C_VCA"]), _VCA_CONCS, "C_VCA")
    r_core_min = float(inputs["r_core_min"])
    r_cutoff = float(inputs["r_cutoff"])

    csv_dir = Path(os.environ["BOXGRAPH_SOURCE"]) / "PLOTCSVFILES"
    if dataset == "arp23":
        prefix = f"Density_Arp23_concsweep_stericR_{a_idx}_r_"
    elif dataset == "vca":
        rval = a_idx * 5 + v_idx
        prefix = f"Density_VCA_Arp23_concsweep_steric_kunbind_0R_{rval}_r_"
    else:
        raise ValueError(f"Unknown dataset setting: {dataset}")

    rep_files = sorted(csv_dir.glob(f"{prefix}*.csv"))
    if not rep_files:
        raise FileNotFoundError(f"No density CSVs matching {prefix}*.csv in {csv_dir}")

    rho_profiles = []
    rae_list = []
    final_times = []

    for fp in rep_files:
        with fp.open(encoding="utf-8") as fh:
            rows = list(csv.reader(fh))
        header = rows[0]
        # header = ['', 'Radius', '0', '1', ..., '<last_t>']
        final_times.append(float(header[-1]))
        radvec = np.array([0.0] + [float(r[1]) for r in rows[1:]], dtype=float)
        dR = radvec[1] - radvec[0]
        Rmax = radvec[-1]
        n_cols = len(header)
        # Exact Cell 39 / Cell 43 tail loop: last n_tail frames of THIS file
        for col_idx in range(n_cols - n_tail, n_cols):
            nvec = np.array([0.0] + [float(r[col_idx]) for r in rows[1:]], dtype=float)
            nmax = float(np.sum(nvec))
            bulkfactor = 3.0 * nmax / (Rmax ** 3)
            with np.errstate(divide="ignore", invalid="ignore"):
                localfactor = (dR ** -1) * (radvec ** -2)
                density_ratio = np.where(radvec > 0.0, nvec * localfactor / bulkfactor, 0.0)
            rho_profiles.append(density_ratio)

            # Exact Cell 39 / Cell 43 R.A.E. (argmax of raw count N, then divide by r^2)
            core_mask = (radvec >= r_core_min) & (radvec < r_cutoff)
            peri_mask = radvec >= r_cutoff
            rc, nc = radvec[core_mask], nvec[core_mask]
            rp, np_ = radvec[peri_mask], nvec[peri_mask]
            ic = int(np.argmax(nc))
            ip = int(np.argmax(np_))
            core_density = nc[ic] / (rc[ic] ** 2)
            peri_density = np_[ip] / (rp[ip] ** 2)
            rae_list.append(float(core_density / peri_density))

    mean_profile = np.mean(rho_profiles, axis=0)
    mask_plot = radvec >= 0.1
    return {
        "rae_mean": float(np.mean(rae_list)),
        "rae_sd": float(np.std(rae_list)),
        "rho_ratio_peak": float(np.max(mean_profile[mask_plot])),
        "n_reps_found": float(len(rep_files)),
        "min_final_time": float(min(final_times)),
    }
