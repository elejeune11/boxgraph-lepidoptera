"""Live software adapter reading PLOTCSVFILES/Fillength_Arp23_concsweep_steric and running Cell 35 converttocdf()."""

import ast
import csv
import json
import math
import os
import types
from pathlib import Path
import numpy as np
import scipy.constants

_ARP_CONCS = [0.0, 10.0, 25.0, 50.0, 75.0, 100.0, 125.0, 150.0]
_CDF_FN = None


def _load_cdf_fn():
    global _CDF_FN
    if _CDF_FN is not None:
        return _CDF_FN
    nb_path = Path(os.environ["BOXGRAPH_SOURCE"]) / "CytoSimanalyze_VASP_Arp_VCA.ipynb"
    nb = json.loads(nb_path.read_text(encoding="utf-8"))
    src = "".join(nb["cells"][35]["source"])
    tree = ast.parse(src)
    fn_node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "converttocdf")
    mod = ast.Module(body=[fn_node], type_ignores=[])
    # NumPy 2.x removed np.trapz in favor of np.trapezoid; provide thin shim without altering notebook code
    np_shim = types.SimpleNamespace(zeros=np.zeros, trapz=getattr(np, "trapz", np.trapezoid))
    ns = {"np": np_shim}
    exec(compile(mod, "CytoSimanalyze_VASP_Arp_VCA.ipynb#cell35", "exec"), ns)
    _CDF_FN = ns["converttocdf"]
    return _CDF_FN


def _arp_index(c_arp):
    for idx, val in enumerate(_ARP_CONCS):
        if abs(c_arp - val) < 1e-6:
            return idx
    raise ValueError(f"C_Arp={c_arp} nM not in Arp23_concsweep_steric grid {_ARP_CONCS}")


def run(inputs, settings):
    cdf_fn = _load_cdf_fn()
    r_idx = _arp_index(float(inputs["C_Arp"]))
    R = float(inputs["R_drop"])
    l_mono = float(inputs["l_mono"])

    csv_dir = Path(os.environ["BOXGRAPH_SOURCE"]) / "PLOTCSVFILES"
    fp = csv_dir / f"Fillength_Arp23_concsweep_stericR_{r_idx}.csv"
    with fp.open(encoding="utf-8") as fh:
        rows = {row[0]: [float(x) for x in row[1:]] for row in csv.reader(fh) if row and row[0]}

    nfil_0 = rows["Nfil_mean"][0]
    nfil_600 = rows["Nfil_mean"][-1]
    lfil_600 = rows["Lfil_mean"][-1]
    lsum_600 = rows["Lsum_mean"][-1]

    vol_L = (4.0 / 3.0) * math.pi * (R ** 3) * 1e-15
    uM_per_um = (1.0 / l_mono) / (scipy.constants.Avogadro * vol_L) * 1e6
    factin_600 = lsum_600 * uM_per_um

    # Execute exact Cell 35 right-edge X grid and converttocdf(X, df[595])
    fp_dist = csv_dir / f"Fillength_dist_Arp23_concsweep_stericR_{r_idx}.csv"
    with fp_dist.open(encoding="utf-8") as fh:
        dist_rows = {row[0]: row[1:] for row in csv.reader(fh) if row and row[0]}
    px_595 = np.array([float(x) for x in dist_rows["595"]], dtype=float)
    X_edges = np.linspace(0.0, 6.5, 100)
    X = X_edges[:-1] + (X_edges[1] - X_edges[0])
    cdfvec = cdf_fn(X, px_595)
    minpos = int(np.argmin(np.abs(cdfvec - 0.5)))

    return {
        "N_fil_600": float(nfil_600),
        "N_fil_fold": float(nfil_600 / nfil_0),
        "L_fil_mean_600": float(lfil_600),
        "L_sum_600": float(lsum_600),
        "F_actin_conc_600": float(factin_600),
        "L_fil_median_cdf": float(X[minpos]),
    }
