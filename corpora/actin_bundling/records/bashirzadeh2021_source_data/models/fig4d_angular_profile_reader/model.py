"""Live reader auditing confocal angular intensity profiles in Fig4c-d.mat (Fig. 4d)."""
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
        1: ("ValuesAlphaoval", "ValuesFascinoval", "ValuesActinoval"),     # Fig. 4d primary circle
        2: ("ValuesAlpha2oval", "ValuesFascin2oval", "ValuesActin2oval"),  # Fig. 4d secondary circle
    }
    k_a, k_f, k_c = key_map[pidx]
    arr_a = np.asarray(mat[k_a], dtype=float)
    arr_f = np.asarray(mat[k_f], dtype=float)
    arr_c = np.asarray(mat[k_c], dtype=float)

    ya = arr_a[:, 1]
    yf = arr_f[:, 1]
    yc = arr_c[:, 1]

    return {
        "n_points": float(len(arr_a)),
        "angular_span": float(arr_a[-1, 0]),
        "alpha_mean_raw": float(ya.mean()),
        "fascin_mean_raw": float(yf.mean()),
        "actin_mean_raw": float(yc.mean()),
        "alpha_cv": float(ya.std(ddof=1) / ya.mean()),
        "fascin_cv": float(yf.std(ddof=1) / yf.mean()),
        "actin_cv": float(yc.std(ddof=1) / yc.mean()),
        "corr_alpha_fascin": float(np.corrcoef(ya, yf)[0, 1]),
        "corr_actin_alpha": float(np.corrcoef(yc, ya)[0, 1]),
        "corr_actin_fascin": float(np.corrcoef(yc, yf)[0, 1]),
    }
