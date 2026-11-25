import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve


def _precompute_phi_table(A, K, zeta, lam, d, R_min, R_max, n_R):
    alpha = (zeta - 2.0 * lam) / K
    nodes, weights = np.polynomial.legendre.leggauss(48)

    def psi_g(phi):
        half = 0.5 * phi
        xs_q = half * (nodes + 1.0)
        ex = np.exp(alpha * xs_q)
        return float(half * np.dot(weights, ex)), float(half * np.dot(weights, A * (xs_q**3 - xs_q) * ex))

    def res_flat(x):
        p1, p2 = float(x[0]), float(x[1])
        mu1, mu2 = A * (p1**3 - p1), A * (p2**3 - p2)
        psi1, g1 = psi_g(p1)
        psi2, g2 = psi_g(p2)
        return [mu1 - mu2, (psi1 * mu1 - g1) - (psi2 * mu2 - g2)]

    p1, p2 = fsolve(res_flat, [-1.0, 1.0], xtol=1e-12)
    p1, p2 = float(p1), float(p2)
    if p1 > p2:
        p1, p2 = p2, p1
    mu_bar = A * (p1**3 - p1)
    psi1, g1 = psi_g(p1)
    P1 = psi1 * mu_bar - g1

    n_grid = 200
    xs0 = np.linspace(p1, p2, n_grid)
    w0 = np.empty_like(xs0)
    for i, x in enumerate(xs0):
        psix, gx = psi_g(float(x))
        w0[i] = max(0.0, -(2.0 / K) * ((psix * mu_bar - gx) - P1) * math.exp(-alpha * float(x)))
    sw0 = np.sqrt(w0)
    I0 = float(np.trapezoid(sw0, xs0))
    S0 = math.exp(alpha * p2) * I0
    S1 = float(np.trapezoid(sw0 * np.exp(alpha * xs0), xs0))
    if abs(zeta - 2.0 * lam) < 1e-8:
        I_phi = float(np.trapezoid(xs0 * sw0, xs0))
        sigma0 = K * I0 + zeta * (p2 * I0 - I_phi)
    else:
        sigma0 = (K / (zeta - 2.0 * lam)) * (zeta * S0 - 2.0 * lam * S1)

    Rs = np.linspace(R_min, R_max, n_R)
    pms = np.empty(n_R)
    pps = np.empty(n_R)
    pm_guess, pp_guess = p1, p2
    for idx in range(n_R - 1, -1, -1):
        R = float(Rs[idx])
        Delta1 = (d - 1.0) * zeta * I0 / R

        def solve_pair(Delta, sigma, init):
            def res(x):
                pm, pp = float(x[0]), float(x[1])
                mu_m = A * (pm**3 - pm)
                mu_p = A * (pp**3 - pp) + Delta
                psim, gm = psi_g(pm)
                psip, gp = psi_g(pp)
                return [mu_m - mu_p, (psip * mu_m - gp) - (psim * mu_m - gm) - (d - 1.0) * sigma / R]
            pm, pp = fsolve(res, init, xtol=1e-11)
            return float(pm), float(pp)

        pm1, pp1 = solve_pair(Delta1, sigma0, [pm_guess, pp_guess])
        muI1 = A * (pm1**3 - pm1)
        xs2 = np.linspace(pm1, pp1, n_grid)
        sw_old = np.interp(xs2, xs0, sw0, left=0.0, right=0.0)
        dx = xs2[1] - xs2[0]
        cum_sw = np.zeros_like(xs2)
        cum_sw[1:] = np.cumsum(0.5 * (sw_old[:-1] + sw_old[1:]) * dx)
        s_vals = (
            (2.0 / K) * (A * (xs2**3 - xs2) - muI1)
            + (2.0 * (d - 1.0) / R) * sw_old
            + (2.0 * (d - 1.0) * zeta / (K * R)) * cum_sw
        )
        integrand = np.exp(alpha * xs2) * s_vals
        cum_int = np.zeros_like(xs2)
        cum_int[1:] = np.cumsum(0.5 * (integrand[:-1] + integrand[1:]) * dx)
        w2 = np.maximum(0.0, np.exp(-alpha * xs2) * cum_int)
        sw2 = np.sqrt(w2)
        I0_2 = float(np.trapezoid(sw2, xs2))
        S0_2 = math.exp(alpha * pp1) * I0_2
        S1_2 = float(np.trapezoid(sw2 * np.exp(alpha * xs2), xs2))
        Delta2 = (d - 1.0) * zeta * I0_2 / R
        if abs(zeta - 2.0 * lam) < 1e-8:
            I_phi2 = float(np.trapezoid(xs2 * sw2, xs2))
            sigma2 = K * I0_2 + zeta * (pp1 * I0_2 - I_phi2)
        else:
            sigma2 = (K / (zeta - 2.0 * lam)) * (zeta * S0_2 - 2.0 * lam * S1_2)
        pm2, pp2 = solve_pair(Delta2, sigma2, [pm1, pp1])
        pms[idx], pps[idx] = pm2, pp2
        pm_guess, pp_guess = pm2, pp2
    return Rs, pms, pps


def run(inputs, settings):
    A = float(inputs["A"])
    K = float(inputs["K"])
    zeta = float(inputs["zeta"])
    lam = float(inputs["lam"])
    d = float(inputs["d"])
    R1_0 = float(inputs["R1_0"])
    R2_0 = float(inputs["R2_0"])
    t_end = float(inputs["t_end"])

    n_R = int(settings.get("n_R_table", 40))
    rtol = float(settings.get("rtol", 1e-6))
    atol = float(settings.get("atol", 1e-8))

    R_min = 3.0
    R_max = max(R1_0, R2_0) * 1.5 + 5.0
    Rs, pms, pps = _precompute_phi_table(A, K, zeta, lam, d, R_min, R_max, n_R)

    def fprime(phi):
        return A * (phi**3 - phi)

    def rhs(t, y):
        R1, R2 = float(y[0]), float(y[1])
        if R1 <= R_min or R2 <= R_min:
            return [0.0, 0.0]
        pm1 = float(np.interp(R1, Rs, pms))
        pp1 = float(np.interp(R1, Rs, pps))
        pm2 = float(np.interp(R2, Rs, pms))
        pp2 = float(np.interp(R2, Rs, pps))
        dmu = fprime(pm2) - fprime(pm1)
        dR1 = dmu / ((pp1 - pm1) * R1)
        dR2 = -dmu / ((pp2 - pm2) * R2)
        return [dR1, dR2]

    def stop_evaporated(t, y):
        return min(y[0], y[1]) - R_min

    stop_evaporated.terminal = True
    stop_evaporated.direction = -1

    sol = solve_ivp(rhs, [0.0, t_end], [R1_0, R2_0], method="RK45", rtol=rtol, atol=atol, events=stop_evaporated)
    R1_f = float(sol.y[0, -1])
    R2_f = float(sol.y[1, -1])
    # For R1_0 > R2_0: forward Ostwald increases (R1 - R2), reverse Ostwald decreases (R1 - R2)
    ostwald_sign = 1.0 if (R1_f - R2_f) > (R1_0 - R2_0) else -1.0
    return {
        "R1_final": R1_f,
        "R2_final": R2_f,
        "dR1": R1_f - R1_0,
        "dR2": R2_f - R2_0,
        "ostwald_sign": ostwald_sign,
    }
