"""Live software adapter executing config-generator expressions from Cells 15, 17, 19 of CytoSimanalyze_VASP_Arp_VCA.ipynb."""

import json
import math
import os
from pathlib import Path
import numpy as np
import scipy.constants

_LINES = None


def _load_lines():
    global _LINES
    if _LINES is not None:
        return _LINES
    nb_path = Path(os.environ["BOXGRAPH_SOURCE"]) / "CytoSimanalyze_VASP_Arp_VCA.ipynb"
    nb = json.loads(nb_path.read_text(encoding="utf-8"))
    c17 = "".join(nb["cells"][17]["source"]).splitlines()
    c19 = "".join(nb["cells"][19]["source"]).splitlines()
    tp_line = next(ln.strip() for ln in c17 if ln.strip().startswith("total_polymer ="))
    gs_line = next(ln.strip() for ln in c17 if ln.strip().startswith("vca_growingspeedvec="))
    ml_line = next(ln.strip() for ln in c17 if ln.strip().startswith("vca_max_lengthvec="))
    sep_line = next(ln.strip() for ln in c19 if "separation =" in ln and "0.9*2*blobsizevec" in ln)
    _LINES = (tp_line, gs_line, ml_line, sep_line)
    return _LINES


def run(inputs, settings):
    tp_line, _, _, _ = _load_lines()
    R = float(inputs["R_drop"])
    nfil = float(inputs["N_fil"])
    l0 = float(inputs["L0_per_fil"])
    kbind = float(inputs["k_bind_vca"])
    nvca = float(inputs["N_vca"])
    l_mono = float(inputs["l_mono"])
    gfac = float(inputs["gmax_factor"])
    r_blob = float(inputs["r_blob"])

    # Execute exact notebook total_polymer assignment
    ns = {"np": np, "Nfil": [nfil], "radiusstring": [R]}
    exec(compile(tp_line, "CytoSimanalyze_VASP_Arp_VCA.ipynb#cell17", "exec"), ns)
    total_polymer = float(ns["total_polymer"][0])

    # Execute exact notebook VCA formulas and separation = 0.9*2*r_blob
    k_bind_fiber = float(np.multiply([kbind], [nvca])[0])
    L_Gmax = float(np.multiply([l_mono * gfac], [nvca])[0])
    L_max_fil = 2.0 * math.pi * R
    L_max_total = nfil * L_max_fil
    site_sep = 0.9 * 2.0 * r_blob

    # Convert polymer length (um) to molar concentration (uM) in spherical volume V = (4/3)*pi*R^3
    vol_L = (4.0 / 3.0) * math.pi * (R ** 3) * 1e-15
    uM_per_um = (1.0 / l_mono) / (scipy.constants.Avogadro * vol_L) * 1e6
    T_actin_conc = total_polymer * uM_per_um
    F_actin_init_conc = (nfil * l0) * uM_per_um

    return {
        "total_polymer": total_polymer,
        "L_max_fil": L_max_fil,
        "L_max_total": L_max_total,
        "T_actin_conc": T_actin_conc,
        "F_actin_init_conc": F_actin_init_conc,
        "k_bind_fiber": k_bind_fiber,
        "L_Gmax": L_Gmax,
        "site_separation": site_sep,
    }
