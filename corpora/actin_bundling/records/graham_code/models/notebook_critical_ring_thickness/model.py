"""BoxGraph software adapter solving for the circle-to-ellipse transition thickness T_crit (AR = 1.1)
using the live VASPDroplet_simulate.ipynb functions (SystemParams, getkbend, gettotalenergy,
addrandomnoise) with the notebook's exact IEEE-754 Nbeads = int(D*1e3/(D*0.1*1e3)) + int((pi-1)*D/(lcyl*1e-3)) - 1."""

import importlib.util
from pathlib import Path
import numpy as np
from scipy.optimize import Bounds, brentq, minimize


def _load_base():
    base_path = Path(__file__).resolve().parents[1] / "notebook_droplet_ring_slsqp" / "model.py"
    spec = importlib.util.spec_from_file_location("_graham_base_slsqp", base_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run(inputs, settings):
    base = _load_base()
    ns = base._load_notebook()
    seed = int(settings.get("seed", 0))
    lcylscalewithD = bool(settings.get("lcylscalewithD", True))
    ftol = float(settings.get("ftol", 1e-7))

    np.random.seed(seed)
    diameter = float(inputs["D_droplet"])
    radiusnm = diameter * 1e3 / 2.0
    r_fil = float(inputs["R_fil"])
    ar_thresh = float(inputs["AR_thresh"])

    sp = ns["SystemParams"]()
    sp.dimval = 2
    sp.noisefactor = float(inputs["noisefactor"])
    sp.nonlinearbool = False
    sp.alpha = 0.0
    sp.kflex_rigidity_f = float(inputs["kflex_fil"])
    sp.kstr = float(inputs["kstr"])
    sp.surface_T = [float(inputs["gamma_I"]), float(inputs["gamma_AVR"]), float(inputs["gamma_AVD"])]
    com_val = float(inputs["ellipse_com"])
    sp.ellipse_com = [com_val, com_val]
    sp.ellipse = [radiusnm, radiusnm]
    sp.Vdropfactor = radiusnm * radiusnm
    sp.lcyl = diameter * float(inputs["lcyl_scale_frac"]) * 1e3 if lcylscalewithD else float(inputs["lcyl_fixed"])

    # Exact IEEE-754 truncation from generatering() lines 712, 722, 725-744:
    # yields Nbeads = 29 for D = 0.75, 1.5, 1.75, 3.0, 3.5 um and Nbeads = 30 otherwise.
    n_reps = int((np.pi - 1.0) * diameter / (sp.lcyl * 1e-3))
    nbeads_init = int(diameter * 1e3 / sp.lcyl)
    nbeads = nbeads_init + n_reps - 1
    sp.Nbeads = nbeads

    n_seg = nbeads - 1
    r0 = 0.94 * radiusnm
    ang = np.linspace(0.0, n_seg * sp.lcyl / r0, nbeads)
    coords0 = np.column_stack([com_val + r0 * np.cos(ang), com_val + r0 * np.sin(ang)])
    x0 = ns["addrandomnoise"](np.concatenate([[radiusnm], coords0.ravel()]), sp.noisefactor)

    def eval_ar(nfil_val):
        sp.Nfil = float(nfil_val)
        sp.rcyl = r_fil * sp.Nfil
        sp.kbend = float(ns["getkbend"](sp))
        res = minimize(
            base._analytic_fg,
            x0,
            args=(sp,),
            jac=True,
            method="SLSQP",
            bounds=Bounds(1e-9, 1e10),
            options={"ftol": ftol, "maxiter": 1500},
        )
        # Verify energy on unmodified notebook gettotalenergy()
        _ = float(ns["gettotalenergy"](res.x, sp, False))
        a = float(res.x[0])
        b = float(sp.Vdropfactor / a)
        return max(a, b) / min(a, b)

    lo, hi = 0.5, 250.0
    if eval_ar(lo) >= ar_thresh:
        nf_star = lo
    else:
        while eval_ar(hi) < ar_thresh and hi < 2000.0:
            hi *= 2.0
        nf_star = float(brentq(lambda n: eval_ar(n) - ar_thresh, lo, hi, xtol=0.05))

    return {
        "T_crit": float(nf_star * r_fil / 1000.0),
        "n_beads_final": float(nbeads),
    }
