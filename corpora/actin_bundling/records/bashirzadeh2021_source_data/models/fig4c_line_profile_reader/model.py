"""Live reader auditing confocal line intensity profiles in Fig4c-d.mat (Fig. 4c & Supp. Fig. 8e)."""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
from scipy.io import loadmat


def run(inputs, settings):
    src = Path(os.environ["BOXGRAPH_SOURCE"])
    mat = loadmat(src / "Fig4c-d.mat")

    pidx = int(round(float(inputs["profile"])))
    key_map = {
        1: ("ValuesAlpha", "ValuesFascin", "ValuesActin"),           # Fig. 4c primary scan
        2: ("ValuesAlpha2", "ValuesFascin2", "ValuesActin2"),        # Fig. 4c secondary scan
        3: ("ValuesAlphaRing", "ValuesFascinRing", "ValuesActinRing"),  # Supp. Fig. 8e actin ring
    }
    k_a, k_f, k_c = key_map[pidx]
    arr_a = np.asarray(mat[k_a], dtype=float)
    arr_f = np.asarray(mat[k_f], dtype=float)
    arr_c = np.asarray(mat[k_c], dtype=float)

    n_af = min(len(arr_a), len(arr_f))
    n_ac = min(len(arr_a), len(arr_c))
    n_fc = min(len(arr_f), len(arr_c))

    ya = arr_a[:, 1]
    yf = arr_f[:, 1]
    yc = arr_c[:, 1]

    return {
        "n_points": float(len(arr_a)),
        "spatial_span": float(arr_a[-1, 0]),
        "alpha_mean_raw": float(ya.mean()),
        "fascin_mean_raw": float(yf.mean()),
        "actin_mean_raw": float(yc.mean()),
        "alpha_norm_peak": float(ya.max() / ya.mean()),
        "fascin_norm_peak": float(yf.max() / yf.mean()),
        "actin_norm_peak": float(yc.max() / yc.mean()),
        "alpha_peak_position": float(arr_a[int(np.argmax(ya)), 0]),
        "fascin_peak_position": float(arr_f[int(np.argmax(yf)), 0]),
        "corr_alpha_fascin": float(np.corrcoef(ya[:n_af], yf[:n_af])[0, 1]),
        "corr_actin_alpha": float(np.corrcoef(yc[:n_ac], ya[:n_ac])[0, 1]),
        "corr_actin_fascin": float(np.corrcoef(yc[:n_fc], yf[:n_fc])[0, 1]),
    }
