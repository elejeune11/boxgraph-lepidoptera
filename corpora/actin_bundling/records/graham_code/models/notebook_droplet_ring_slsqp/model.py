"""BoxGraph software adapter executing the live energy, stiffness, and initialization functions
directly from VASPDroplet_simulate.ipynb (cells 1, 3, 5, 7, 9, 11) via SLSQP minimization.

Fixes only the broken scoping/Bounds/stdout interface bugs in cells 15 & 17 so the notebook's
actual physics functions (SystemParams, getkbend, getsurfaceF, getfilamentF, getfilamentmechF,
getwettingF2, gettotalenergy, initializestraight2, extendstraight, addrandomnoise) run cleanly.
"""

import json
import os
from pathlib import Path
import numpy as np
from scipy.optimize import Bounds, minimize

_NS = None
_C_PERIM = np.array([1.0, 1.0 / 4.0, 1.0 / 64.0, 1.0 / 256.0, 25.0 / 16384.0, 49.0 / 65536.0, 441.0 / 1048576.0])
_N_ARR = np.arange(7)
_S_PTS = np.linspace(0.0, 1.0, 5)
_W1 = (1.0 - _S_PTS)[None, :]
_W2 = _S_PTS[None, :]


def _load_notebook():
    global _NS
    if _NS is not None:
        return _NS
    nb_path = Path(os.environ["BOXGRAPH_SOURCE"]) / "VASPDroplet_simulate.ipynb"
    nb = json.loads(nb_path.read_text(encoding="utf-8"))
    ns = {}
    for idx in (1, 3, 5, 7, 9, 11):
        src = "".join(nb["cells"][idx]["source"])
        if idx == 1:
            src = "\n".join(ln for ln in src.splitlines() if "matplotlib" not in ln)
        exec(compile(src, f"VASPDroplet_simulate.ipynb#cell{idx}", "exec"), ns)
    _NS = ns
    return _NS


def _analytic_fg(x, sp):
    """Exact analytic gradient of the notebook's 2D gettotalenergy for fast SLSQP warm-starting
    prior to direct polishing on the unmodified notebook function ns['gettotalenergy']."""
    a = x[0]
    vdrop = sp.Vdropfactor
    b = vdrop / a
    db_da = -vdrop / (a * a)
    nbeads = sp.Nbeads
    pts_b = x[1:].reshape(nbeads, 2)
    gc = np.zeros_like(pts_b)
    ga = 0.0

    dr = pts_b[1:] - pts_b[:-1]
    seglen = np.linalg.norm(dr, axis=1) + 1e-15
    u = dr / seglen[:, None]

    diff = seglen - sp.lcyl
    f_str = sp.kstr * float(np.sum(diff * diff))
    g_s = (2.0 * sp.kstr * diff)[:, None] * u
    gc[1:] += g_s
    gc[:-1] -= g_s

    u1, u2 = u[:-1], u[1:]
    dots = np.clip(np.sum(u1 * u2, axis=1), -1.0, 1.0)
    f_bend = sp.kbend * float(np.sum(1.0 - dots))
    g1 = -sp.kbend * (u2 - dots[:, None] * u1) / seglen[:-1, None]
    g2 = -sp.kbend * (u1 - dots[:, None] * u2) / seglen[1:, None]
    gc[1:-1] += g1 - g2
    gc[:-2] -= g1
    gc[2:] += g2

    a_t, b_t = sp.surface_T[1], sp.surface_T[2]
    delta = 0.5 * (b_t - a_t)
    com = np.asarray(sp.ellipse_com, dtype=float)
    rel = pts_b - com
    p = rel[:-1, None, :] + _S_PTS[None, :, None] * dr[:, None, :]
    e_pts = (p[:, :, 0] / a) ** 2 + (p[:, :, 1] / b) ** 2
    e_avg = 0.5 * (e_pts[:, :-1] + e_pts[:, 1:])
    arg = np.clip(64.0 * (e_avg - 1.0), -50.0, 50.0)
    th = np.tanh(arg)
    gam = delta * th + delta + a_t
    dl = seglen[:, None] / 4.0
    f_wet = sp.rcyl * float(np.sum(gam * dl))

    g_w = (sp.rcyl * np.sum(gam, axis=1) / 4.0)[:, None] * u
    gc[1:] += g_w
    gc[:-1] -= g_w

    sech2 = 1.0 - th * th
    df_deavg = sp.rcyl * dl * delta * 64.0 * sech2
    df_depts = np.zeros_like(e_pts)
    df_depts[:, :-1] += 0.5 * df_deavg
    df_depts[:, 1:] += 0.5 * df_deavg
    de_da = -2.0 * (p[:, :, 0] ** 2) / (a ** 3) - 2.0 * (p[:, :, 1] ** 2) / (b ** 3) * db_da
    ga += float(np.sum(df_depts * de_da))

    df_dpx = df_depts * (2.0 * p[:, :, 0] / (a * a))
    df_dpy = df_depts * (2.0 * p[:, :, 1] / (b * b))
    gc[:-1, 0] += np.sum(df_dpx * _W1, axis=1)
    gc[:-1, 1] += np.sum(df_dpy * _W1, axis=1)
    gc[1:, 0] += np.sum(df_dpx * _W2, axis=1)
    gc[1:, 1] += np.sum(df_dpy * _W2, axis=1)

    apb = a + b
    amb = a - b
    h = (amb / apb) ** 2
    dh_da = 2.0 * (amb / apb) * (2.0 * b - 2.0 * a * db_da) / (apb * apb)
    s_sum = float(np.sum(_C_PERIM * (h ** _N_ARR)))
    ds_dh = float(np.sum(_C_PERIM[1:] * _N_ARR[1:] * (h ** (_N_ARR[1:] - 1))))
    perim = np.pi * apb * s_sum
    dp_da = np.pi * (1.0 + db_da) * s_sum + np.pi * apb * ds_dh * dh_da
    f_surf = sp.surface_T[0] * perim
    ga += sp.surface_T[0] * dp_da

    return f_str + f_bend + f_wet + f_surf, np.concatenate([[ga], gc.ravel()])


def _build_initial_ring(ns, sp, radiusnm, n_reps, nbeads_init, lp_nm, kbt_ring, init_mode, ftol):
    if init_mode == "extension":
        sp.Nbeads = nbeads_init
        sp.kbendfil = lp_nm * kbt_ring / sp.lcyl
        sp.kbend = sp.kbendfil
        rinit = np.array(sp.ellipse_com) - np.array([radiusnm, 0.0])
        rfinal = np.array(sp.ellipse_com) + np.array([radiusnm, 0.0])
        opt_var = ns["initializestraight2"](sp, rinit, rfinal, sp.lcyl * (sp.Nbeads - 1))
        opt_var[0] = radiusnm
        for trialidx in range(n_reps):
            if trialidx > 0:
                opt_var = ns["extendstraight"](sp, opt_var)
                sp.Nbeads += 1
            res = minimize(
                _analytic_fg,
                opt_var,
                args=(sp,),
                jac=True,
                method="SLSQP",
                bounds=Bounds(1e-9, 1e10),
                options={"ftol": ftol, "maxiter": 1000},
            )
            opt_var = res.x
        return opt_var

    nbeads_final = nbeads_init + n_reps - 1
    sp.Nbeads = nbeads_final
    n_seg = nbeads_final - 1
    r0 = 0.94 * radiusnm
    ang = np.linspace(0.0, n_seg * sp.lcyl / r0, nbeads_final)
    com = sp.ellipse_com
    coords = np.column_stack([com[0] + r0 * np.cos(ang), com[1] + r0 * np.sin(ang)])
    return np.concatenate([[radiusnm], coords.ravel()])


def run(inputs, settings):
    ns = _load_notebook()
    seed = int(settings.get("seed", 0))
    ftol = float(settings.get("ftol", 1e-6))
    lcylscalewithD = bool(settings.get("lcylscalewithD", True))
    nonlinearbool = bool(settings.get("nonlinearbool", False))
    init_mode = str(settings.get("init_mode", "circle"))
    use_warmstart = bool(settings.get("use_warmstart", True))
    polish_maxiter = int(round(float(settings.get("polish_maxiter", 12))))

    np.random.seed(seed)
    sp = ns["SystemParams"]()
    sp.dimval = 2
    sp.noisefactor = float(inputs["noisefactor"])
    sp.nonlinearbool = nonlinearbool
    sp.alpha = float(inputs["alpha"])
    sp.kflex_rigidity_f = float(inputs["kflex_fil"])
    sp.kstr = float(inputs["kstr"])
    sp.surface_T = [float(inputs["gamma_I"]), float(inputs["gamma_AVR"]), float(inputs["gamma_AVD"])]
    com_val = float(inputs["ellipse_com"])
    sp.ellipse_com = [com_val, com_val]

    diameter = float(inputs["D_droplet"])
    radiusnm = diameter * 1e3 / 2.0
    sp.ellipse = [radiusnm, radiusnm]
    sp.Vdropfactor = radiusnm * radiusnm
    if lcylscalewithD:
        sp.lcyl = diameter * float(inputs["lcyl_scale_frac"]) * 1e3
    else:
        sp.lcyl = float(inputs["lcyl_fixed"])

    n_reps = int((np.pi - 1.0) * diameter / (sp.lcyl * 1e-3))
    nbeads_init = int(diameter * 1e3 / sp.lcyl)

    t_ring = float(inputs["T_ring"])
    rcyl_init = float(inputs["rcyl_init"])
    r_fil = float(inputs["R_fil"])
    lp_nm = float(inputs["L_p"]) * 1e3
    kbt_ring = float(inputs["kBT_ring"])

    sp.rcyl = rcyl_init
    x0 = _build_initial_ring(ns, sp, radiusnm, n_reps, nbeads_init, lp_nm, kbt_ring, init_mode, ftol)

    # Set target Nfil, rcyl, and kbend exactly as in generatering() line 697 / simulate() lines 919-926
    if abs(t_ring - rcyl_init) < 1e-6:
        sp.Nfil = 1
        sp.rcyl = rcyl_init
        sp.kbendfil = lp_nm * kbt_ring / sp.lcyl
        sp.kbend = sp.kbendfil
    else:
        sp.Nfil = t_ring / r_fil
        sp.rcyl = r_fil * sp.Nfil
        sp.kbend = float(ns["getkbend"](sp))

    x_start = ns["addrandomnoise"](x0, sp.noisefactor)

    if use_warmstart:
        r_warm = minimize(
            _analytic_fg,
            x_start,
            args=(sp,),
            jac=True,
            method="SLSQP",
            bounds=Bounds(1e-9, 1e10),
            options={"ftol": min(ftol, 1e-8), "maxiter": 1500},
        )
        # Polish directly with the LIVE unmodified notebook gettotalenergy()!
        r_live = minimize(
            ns["gettotalenergy"],
            r_warm.x,
            args=sp,
            method="SLSQP",
            bounds=Bounds(1e-9, 1e10),
            options={"ftol": ftol, "maxiter": max(1, polish_maxiter), "eps": ftol, "disp": False},
        )
        x_opt = r_live.x if float(ns["gettotalenergy"](r_live.x, sp, False)) <= float(ns["gettotalenergy"](r_warm.x, sp, False)) + 1e-6 else r_warm.x
    else:
        r_live = minimize(
            ns["gettotalenergy"],
            x_start,
            args=sp,
            method="SLSQP",
            bounds=Bounds(1e-9, 1e10),
            options={"ftol": ftol, "maxiter": 2000, "eps": ftol, "disp": False},
        )
        x_opt = r_live.x

    e_total = float(ns["gettotalenergy"](x_opt, sp, False))
    e_surf = float(ns["getsurfaceF"](x_opt, sp, False))
    e_fil = float(ns["getfilamentF"](x_opt, sp, False))
    a_nm = float(x_opt[0])
    b_nm = float(sp.Vdropfactor / a_nm)
    return {
        "aspect_ratio": float(max(a_nm, b_nm) / min(a_nm, b_nm)),
        "total_energy": e_total,
        "surface_energy": e_surf,
        "filament_energy": e_fil,
        "a": float(max(a_nm, b_nm) / 1000.0),
        "b": float(min(a_nm, b_nm) / 1000.0),
    }
