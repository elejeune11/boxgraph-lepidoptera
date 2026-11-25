import math
import numpy as np
from scipy.optimize import brentq


def _eval_point(Z, Om, A, Kd, Dc, Dap, Dam, Om0, Omd, Ks, include_uninit_col=True):
    """Exact line-by-line port of 1D/stabilitydiagram.jl:55-126."""
    B = Z  # line 90: Pars["B"] = Pars["Z"]
    Na_HSS = float(
        brentq(
            lambda Na: -(Omd * A / Kd) * Na * Na + Om0 * (1.0 - Na) * (1.0 + Om * Na * Na),
            0.0,
            1.0,
            xtol=1e-14,
            rtol=1e-14,
        )
    )
    Ni_HSS = 1.0 - Na_HSS
    C_HSS = (A / Kd) * Na_HSS
    Nm_HSS = Na_HSS - Ni_HSS
    dPi_HSS = 2.0 * Z * C_HSS - 3.0 * B * (C_HSS ** 2)

    Nk = len(Ks)
    K2 = Ks * Ks
    fac = K2 / (1.0 + K2)

    M = np.empty((Nk, 3, 3), dtype=float)
    M[:, 0, 0] = fac * dPi_HSS * C_HSS - Dc * K2 - Kd
    M[:, 0, 1] = A / 2.0
    M[:, 0, 2] = A / 2.0
    M[:, 1, 0] = fac * dPi_HSS
    M[:, 1, 1] = -(Dap / 2.0) * K2
    M[:, 1, 2] = -(Dam / 2.0) * K2
    M[:, 2, 0] = fac * dPi_HSS * Nm_HSS - (1.0 + Nm_HSS) * Omd
    M[:, 2, 1] = (
        -(Dam / 2.0) * K2
        + Om0
        - 0.25 * (Nm_HSS - 3.0) * (Nm_HSS + 1.0) * Om * Om0
        - C_HSS * Omd
    )
    M[:, 2, 2] = (
        -(Dap / 2.0) * K2
        - 0.25 * (4.0 + (1.0 + Nm_HSS) * (3.0 * Nm_HSS - 1.0) * Om) * Om0
        - C_HSS * Omd
    )

    evals = np.linalg.eigvals(M)  # (Nk, 3)
    if include_uninit_col:
        # Line 52 allocates (3, Nk + 1) zeros while lines 103-110 fill only p = 1:Nk
        evals_all = np.vstack([evals, np.zeros((1, 3), dtype=complex)])
    else:
        evals_all = evals

    idx = np.unravel_index(int(np.argmax(evals_all.real)), evals_all.shape)
    max_re = float(evals_all[idx].real)
    max_im = float(abs(evals_all[idx].imag))
    k_dom = float(Ks[idx[0]]) if idx[0] < Nk else 0.0

    if max_re > 0.0:
        code = 1.0 if max_im < 1e-14 else 2.0
    else:
        code = 0.0
    return Na_HSS, C_HSS, max_re, max_im, k_dom, code


def run(inputs, settings):
    A = float(inputs["A"])
    Kd = float(inputs["Kd"])
    Dc = float(inputs["Dc"])
    Dap = float(inputs["Dap"])
    Dam = float(inputs["Dam"])
    Om0 = float(inputs["Omega0"])
    Omd = float(inputs["Omega_d"])
    Z = float(inputs["Z"])
    Om = float(inputs["Omega"])
    L = float(inputs["L"])
    Nx = int(round(float(inputs["Nx"])))

    include_uninit_col = bool(settings.get("include_uninit_col", True))
    om_max = float(settings.get("omega_max_bisect", 40.0))
    n_iter = int(round(settings.get("bisect_iters", 50)))

    dx = L / Nx
    Kmax = math.pi / dx
    Kmin = 2.0 * math.pi / L
    Nk = int(math.floor((Kmax - Kmin) / Kmin))
    Ks = np.arange(1, Nk + 1, dtype=float) * Kmin

    Na_HSS, C_HSS, max_re, max_im, k_dom, code = _eval_point(
        Z, Om, A, Kd, Dc, Dap, Dam, Om0, Omd, Ks, include_uninit_col=include_uninit_col
    )

    # Also locate exact discrete-K onset Omega_crit(Z) on K_p = p * Kmin
    _, _, r0, _, k0, c0 = _eval_point(Z, 0.0, A, Kd, Dc, Dap, Dam, Om0, Omd, Ks, include_uninit_col=False)
    if c0 > 0.0:
        om_c, k_c = 0.0, k0
    else:
        lo, hi = 0.0, om_max
        for _ in range(n_iter):
            mid = 0.5 * (lo + hi)
            _, _, rm, _, _, _ = _eval_point(Z, mid, A, Kd, Dc, Dap, Dam, Om0, Omd, Ks, include_uninit_col=False)
            if rm <= 0.0:
                lo = mid
            else:
                hi = mid
        _, _, _, _, k_c, _ = _eval_point(Z, hi, A, Kd, Dc, Dap, Dam, Om0, Omd, Ks, include_uninit_col=False)
        om_c = hi

    return {
        "Na_HSS": Na_HSS,
        "C_HSS": C_HSS,
        "MaxEigenvRe": max_re,
        "MaxEigenvIm": max_im,
        "K_dominant": k_dom,
        "regime_code": code,
        "Omega_crit_discrete": om_c,
        "K_crit_discrete": k_c,
    }
