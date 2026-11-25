"""Live software adapter executing interpolateallmonomers(), getradialdensitysnap(), and Cell 27/39 radial shell formulas."""

import ast
import json
import os
from pathlib import Path
import numpy as np

_FNS = None


def _load_fns():
    global _FNS
    if _FNS is not None:
        return _FNS
    nb_path = Path(os.environ["BOXGRAPH_SOURCE"]) / "CytoSimanalyze_VASP_Arp_VCA.ipynb"
    nb = json.loads(nb_path.read_text(encoding="utf-8"))
    src = "".join(nb["cells"][27]["source"])
    tree = ast.parse(src)
    wanted = [
        n for n in tree.body
        if isinstance(n, ast.FunctionDef) and n.name in ("interpolateallmonomers", "getradialdensitysnap")
    ]
    mod = ast.Module(body=wanted, type_ignores=[])
    ns = {"np": np}
    exec(compile(mod, "CytoSimanalyze_VASP_Arp_VCA.ipynb#cell27", "exec"), ns)
    _FNS = (ns["interpolateallmonomers"], ns["getradialdensitysnap"])
    return _FNS


def run(inputs, settings):
    interp_fn, hist_fn = _load_fns()
    L_fil = float(inputs["L_fil"])
    L_seg = float(inputs["L_seg"])
    l_mono = float(inputs["l_mono"])
    R_drop = float(inputs["R_drop"])

    # 1. Live test of interpolateallmonomers() on a straight segmented filament (L_seg = 0.1 um)
    n_seg = max(1, int(round(L_fil / L_seg)))
    xs = np.linspace(0.0, L_fil, n_seg + 1)
    fc = np.column_stack([xs, np.zeros_like(xs), np.zeros_like(xs)])
    pts = interp_fn(fc)
    n_pts = float(pts.shape[0])
    exact_monomers = L_fil / l_mono

    # Also exercise getradialdensitysnap() on a short segment inside the sphere
    edges = np.linspace(0.0, R_drop, 41)
    _ = hist_fn([fc[:2]], edges, np.zeros(len(edges) - 1))

    # 2. Exact Cell 27 Radialbincenter (right bin edge!) + Cell 39 density_ratio on a homogeneous sphere
    dR = edges[1] - edges[0]
    rad_right = edges[:-1] + dR  # exact Cell 27 line 1728 formula
    # Exact analytical counts in each shell for a uniform sphere of total count N_total = 1.0:
    # P(r in [r_lo, r_hi]) = (r_hi^3 - r_lo^3) / R^3
    N_shell = (edges[1:] ** 3 - edges[:-1] ** 3) / (R_drop ** 3)
    N_max = float(np.sum(N_shell))
    bulkfactor = 3.0 * N_max / (R_drop ** 3)  # exact Cell 39 line 2455 formula
    localfactor = (dR ** -1) * (rad_right ** -2)  # exact Cell 39 line 2451 formula
    density_ratio = N_shell * localfactor / bulkfactor

    idx_r01 = int(np.argmin(np.abs(rad_right - 0.1)))
    idx_mid = int(np.argmin(np.abs(rad_right - 0.5)))
    idx_out = len(rad_right) - 1

    return {
        "interpolated_monomers": n_pts,
        "monomer_count_ratio": n_pts / exact_monomers,
        "uniform_rho_ratio_outer": float(density_ratio[idx_out]),
        "uniform_rho_ratio_mid": float(density_ratio[idx_mid]),
        "uniform_rho_ratio_r01": float(density_ratio[idx_r01]),
    }
