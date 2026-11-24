import math
import numpy as np


def run(inputs, settings):
    # Work in nm so k_b (pN*nm^2) and gamma_I (pN) yield energies directly in pN*nm
    R_nm = inputs["R"] * 1000.0
    L_nm = inputs["L"] * 1000.0
    k_b = inputs["k_b"]
    gamma_I = inputs["gamma_I"]
    AR = inputs["AR"]
    n_arc = max(200, int(round(settings.get("n_arc", 4000))))

    a = R_nm * math.sqrt(AR)
    b = R_nm / math.sqrt(AR)

    # Discretize [0, 4*pi] so any arc of length L <= P wraps smoothly around the ellipse
    ts = np.linspace(0.0, 4.0 * math.pi, 2 * n_arc + 1)
    sin_t = np.sin(ts)
    cos_t = np.cos(ts)
    denom = (a * sin_t) ** 2 + (b * cos_t) ** 2
    ds_dt = np.sqrt(denom)
    # kappa(t) = a*b / denom^(3/2), so kappa^2 * ds/dt = (a*b)^2 / denom^(5/2)
    k2_ds_dt = ((a * b) ** 2) / (denom ** 2.5)

    dt = ts[1] - ts[0]
    cum_s = np.concatenate([[0.0], np.cumsum(0.5 * (ds_dt[:-1] + ds_dt[1:]) * dt)])
    cum_k2 = np.concatenate([[0.0], np.cumsum(0.5 * (k2_ds_dt[:-1] + k2_ds_dt[1:]) * dt)])

    P = float(cum_s[n_arc])
    E_surf = gamma_I * P

    # Search starting arc length s0 in [0, P) for minimum bending energy placement
    n_starts = max(200, n_arc // 2)
    s0 = np.linspace(0.0, P, n_starts, endpoint=False)
    k2_start = np.interp(s0, cum_s, cum_k2)
    k2_end = np.interp(s0 + L_nm, cum_s, cum_k2)
    E_bend = float(0.5 * k_b * np.min(k2_end - k2_start))

    return {
        "E_bend": E_bend,
        "E_surf": float(E_surf),
        "P_over_2L": float(P / (2.0 * L_nm)),
    }
