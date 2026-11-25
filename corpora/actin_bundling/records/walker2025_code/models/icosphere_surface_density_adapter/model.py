"""Live BoxGraph software adapter calling icosphere.py surface-occupancy functions directly."""

import importlib.util
import math
import os
from pathlib import Path

import numpy as np

_ICO = None


def _load_icosphere():
    global _ICO
    if _ICO is not None:
        return _ICO
    src = Path(os.environ["BOXGRAPH_SOURCE"]) / "Analysis_Scripts" / "icosphere.py"
    spec = importlib.util.spec_from_file_location("walker2025_icosphere_surf", src)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    _ICO = mod
    return _ICO


def run(inputs, settings):
    ic = _load_icosphere()
    R_um = float(inputs["R_drop"])
    subdiv = int(round(float(inputs["n_subdiv"])))
    r_fil = float(inputs["r_fil"])
    L_fil = float(inputs["L_fil"])
    L_seg = float(inputs["L_seg"])

    meshvar = ic.icosphere(subdiv, R_um)
    n_triangles = int(meshvar[1].shape[0])
    unitnormal = ic.generateunitnormals(meshvar, R_um)

    # Construct an equatorial arc of radius r_fil discretized into segments of chord length L_seg
    n_segs = max(1, int(round(L_fil / L_seg)))
    dtheta = 2.0 * math.asin(min(0.999999, L_seg / (2.0 * r_fil)))
    angles = np.arange(n_segs + 1, dtype=float) * dtheta
    fc = np.column_stack([r_fil * np.cos(angles), r_fil * np.sin(angles), np.zeros_like(angles)])

    interp = ic.interpolateallmonomers(fc)
    actincounter = ic.generatedensityfield(meshvar, unitnormal, [fc], ic.SearchAlgoType.ORIGINAL)

    is_all_zero = np.argwhere(actincounter == 0.0)
    frac_occupied = 1.0 - len(is_all_zero) / float(n_triangles)
    occupied = n_triangles - len(is_all_zero)

    # Effective radial shell thickness hardcoded by dthreshold_4m_center = 0.9*0.9 in icosphere.py:224
    shell_thickness_nm = (R_um - 0.9) * 1000.0

    return {
        "shell_thickness": float(shell_thickness_nm),
        "interpolated_monomers": float(interp.shape[0]),
        "occupied_triangles": float(occupied),
        "frac_occupied": float(frac_occupied),
    }
