import math
import numpy as np
from scipy.optimize import brentq


def run(inputs, settings):
    D_B = inputs["D_B"]
    gb = inputs["gamma_beta"]
    k_f = inputs["k_f"]
    k_b = inputs["k_b"]
    V_s = inputs["V_s"]
    phi_B_bar = inputs["phi_B_bar"]
    N = int(round(inputs["N"]))
    rtol = float(settings.get("rtol", 1.0e-6))

    passive = (k_f == 0.0 and k_b == 0.0)

    def phi0B_of(R, m):
        if passive:
            Vtot = m * (4.0 * math.pi / 3.0) * R**3
            if Vtot >= V_s:
                return -1.0
            return (V_s * phi_B_bar - Vtot) / (V_s - Vtot)
        return (k_f + 4.0 * math.pi * m * D_B * gb / V_s) / (
            k_f + k_b + 4.0 * math.pi * m * D_B * R / V_s
        )

    def find_roots_m(m):
        if passive:
            R_hi = (3.0 * V_s * phi_B_bar / (4.0 * math.pi * m)) ** (1.0 / 3.0) * 0.999999
        else:
            R_hi = max(
                math.sqrt(3.0 * D_B / k_b) * 2.0,
                (3.0 * V_s / (4.0 * math.pi * m)) ** (1.0 / 3.0),
            )

        def f(R):
            return D_B * (phi0B_of(R, m) / R - gb / R**2) - k_b * R / 3.0

        R_lo = max(gb * 0.01, R_hi * 1e-6)
        grid = np.geomspace(R_lo, R_hi, 500)
        vals = [f(float(r)) for r in grid]
        roots = []
        for i in range(len(grid) - 1):
            if vals[i] * vals[i + 1] <= 0.0:
                if vals[i] == 0.0:
                    roots.append(float(grid[i]))
                elif vals[i + 1] != 0.0:
                    roots.append(
                        float(brentq(f, float(grid[i]), float(grid[i + 1]), xtol=rtol * R_lo, rtol=rtol))
                    )
        return sorted(roots)

    def exact_jacobian_evals(R_bar, m):
        p0 = phi0B_of(R_bar, m)
        # Explicit (m+1)x(m+1) Jacobian matrix of Eqs. (14) and (13) with X = [phi_0^B, R_1, ..., R_m]
        J = np.zeros((m + 1, m + 1), dtype=float)
        J[0, 0] = -(k_f + k_b) - 4.0 * math.pi * D_B * m * R_bar / V_s
        Jii = D_B * (2.0 * gb / R_bar**3 - p0 / R_bar**2) - k_b / 3.0
        for i in range(1, m + 1):
            J[0, i] = -4.0 * math.pi * D_B * p0 / V_s
            J[i, 0] = D_B / R_bar
            J[i, i] = Jii
        evals = np.linalg.eigvals(J)
        # Collective 2x2 block eigenvalues (exact same non-degenerate pair as full Jacobian)
        M2 = np.array([[J[0, 0], m * J[0, 1]], [J[1, 0], Jii]], dtype=float)
        ev2 = np.linalg.eigvals(M2)
        if np.max(np.abs(ev2.imag)) > 1e-12 * max(1.0, float(np.max(np.abs(ev2.real)))):
            e_plus_real = None
        else:
            e_plus_real = float(np.max(ev2.real))
        return float(Jii), e_plus_real, float(np.max(evals.real)), p0

    # Evaluate requested N-droplet symmetric state
    rN = find_roots_m(N)
    if len(rN) >= 2:
        R_sym = rN[-1]
        lam_sym, e_plus_exact, _, phi_0_B_sym = exact_jacobian_evals(R_sym, N)
    else:
        R_sym = None
        phi_0_B_sym = None
        lam_sym = None
        e_plus_exact = None

    # Pair (N=2) phase portrait fixed-point classification (Fig. 4)
    r1 = find_roots_m(1)
    r2 = find_roots_m(2)
    n_stable_pair = 1  # (0, 0) is stable whenever R_crit > 0
    n_small_unstable_pair = 0
    n_large_sym_unstable_pair = 0
    if len(r1) == 2:
        n_small_unstable_pair += 2  # (R_crit1, 0) and (0, R_crit1)
        n_stable_pair += 2  # (R_bar1, 0) and (0, R_bar1)
    if len(r2) == 2:
        n_small_unstable_pair += 1  # (R_crit2, R_crit2)
        Jii2, _, maxev2, _ = exact_jacobian_evals(r2[1], 2)
        if passive:
            n_large_sym_unstable_pair = 1
        elif maxev2 < 0.0 and Jii2 < 0.0:
            n_stable_pair += 1
        else:
            n_large_sym_unstable_pair = 1

    # Maximal stable droplet count N_max in finite system (Fig. 6a)
    N_max = 0
    if not passive:
        for m in range(1, 500):
            rm = find_roots_m(m)
            if len(rm) < 2:
                break
            Jii_m, _, maxev_m, _ = exact_jacobian_evals(rm[1], m)
            if maxev_m < 0.0 and Jii_m < 0.0:
                N_max = m
            else:
                break
    else:
        N_max = 1 if len(r1) == 2 else 0

    return {
        "R_sym": R_sym,
        "phi_0_B_sym": phi_0_B_sym,
        "lam_sym": lam_sym,
        "e_plus_exact": e_plus_exact,
        "n_stable_pair": float(n_stable_pair),
        "n_small_unstable_pair": float(n_small_unstable_pair),
        "n_large_sym_unstable_pair": float(n_large_sym_unstable_pair),
        "N_max": float(N_max),
    }
