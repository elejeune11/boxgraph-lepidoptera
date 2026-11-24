import math
import numpy as np
from scipy.optimize import brentq, root


def run(inputs, settings):
    D_B = inputs["D_B"]
    gb = inputs["gamma_beta"]
    psi = inputs["psi"]
    k_f = inputs["k_f"]
    k_b = inputs["k_b"]
    k = inputs["k"]
    phibar = inputs["phibar"]
    V_s = inputs["V_s"]
    Q_1 = inputs["Q_1"]
    Q_2 = inputs["Q_2"]
    rtol = float(settings.get("rtol", 1.0e-6))

    k_s = k_f + k * psi

    def bg_ss(R1, R2):
        Vtot = (4.0 * math.pi / 3.0) * (R1**3 + R2**3)
        if Vtot >= 0.999999 * V_s:
            return None, None
        A0 = phibar - (Vtot / V_s) * psi
        cB = (V_s - Vtot) / V_s
        n_act = (1 if R1 > 0.0 else 0) + (1 if R2 > 0.0 else 0)
        R_sum = (R1 if R1 > 0.0 else 0.0) + (R2 if R2 > 0.0 else 0.0)
        qa = -k * cB
        qb = -k_f * cB - k_b + k * A0 - 4.0 * math.pi * D_B * R_sum / (V_s - Vtot)
        qc = k_f * A0 + 4.0 * math.pi * D_B * n_act * gb / (V_s - Vtot)
        if abs(qa) < 1e-30:
            pB = -qc / qb if qb != 0.0 else 0.0
        else:
            disc = qb**2 - 4.0 * qa * qc
            if disc < 0.0:
                return None, None
            pB = (-qb - math.sqrt(disc)) / (2.0 * qa)
            if pB < 0.0:
                pB = 0.0
        pA = A0 - cB * pB
        return pA, pB

    def dV_pair(R1, R2):
        pA, pB = bg_ss(R1, R2)
        if pA is None:
            return None, None, None, None
        V1 = (4.0 * math.pi / 3.0) * R1**3
        V2 = (4.0 * math.pi / 3.0) * R2**3
        dV1 = (4.0 * math.pi * D_B / psi) * (pB * R1 - gb) + V1 * (k_s * pA / psi - k_b) + Q_1 * pA / psi
        dV2 = (4.0 * math.pi * D_B / psi) * (pB * R2 - gb) + V2 * (k_s * pA / psi - k_b) + Q_2 * pA / psi
        return dV1, dV2, pA, pB

    R_scale = (3.0 * V_s * max(phibar, 1e-6) / (4.0 * math.pi * psi)) ** (1.0 / 3.0)
    R_lo = max(1e-4 * gb, 1e-6 * R_scale)
    R_hi = min(1.2 * R_scale, (3.0 * 0.99 * V_s / (4.0 * math.pi)) ** (1.0 / 3.0))
    grid = np.geomspace(R_lo, R_hi, 500)

    # Single-droplet steady states with core Q_1 (R2 = 0)
    vals1 = [dV_pair(float(r), 0.0)[0] for r in grid]
    roots1 = []
    for i in range(len(grid) - 1):
        if vals1[i] is not None and vals1[i + 1] is not None and vals1[i] * vals1[i + 1] < 0.0:
            roots1.append(
                float(brentq(lambda r: dV_pair(r, 0.0)[0], float(grid[i]), float(grid[i + 1]), xtol=rtol * R_lo, rtol=rtol))
            )
    R_single_bar = roots1[-1] if roots1 else None

    q_scale = max(abs(Q_1), abs(Q_2))
    is_symmetric = (q_scale == 0.0) or (abs(Q_1 - Q_2) <= 1e-12 * q_scale)

    if is_symmetric:
        vals2 = [dV_pair(float(r), float(r))[0] for r in grid]
        roots2 = []
        for i in range(len(grid) - 1):
            if vals2[i] is not None and vals2[i + 1] is not None and vals2[i] * vals2[i + 1] < 0.0:
                roots2.append(
                    float(brentq(lambda r: dV_pair(r, r)[0], float(grid[i]), float(grid[i + 1]), xtol=rtol * R_lo, rtol=rtol))
                )
        R1_pair = roots2[-1] if roots2 else None
        R2_pair = R1_pair
        V1_over_V2 = 1.0 if R1_pair is not None else None
        lam_pair = None
        stable2 = False
        if R1_pair is not None:
            h = R1_pair * 1e-5
            dV1p, dV2p, _, _ = dV_pair(R1_pair + h, R1_pair - h)
            dR1p = dV1p / (4.0 * math.pi * (R1_pair + h) ** 2)
            dR2p = dV2p / (4.0 * math.pi * (R1_pair - h) ** 2)
            lam_pair = float((dR1p - dR2p) / (2.0 * h))
            dV1s, _, _, _ = dV_pair(R1_pair + h, R1_pair + h)
            lam_tot = float(dV1s / (4.0 * math.pi * (R1_pair + h) ** 2 * h))
            stable2 = (lam_pair < 0.0 and lam_tot < 0.0)

        if Q_1 == 0.0 and Q_2 == 0.0:
            n_stable_states = 1 + (2 if len(roots1) == 2 else 0) + (1 if stable2 else 0)
            n_crit_radii = (2 if len(roots1) == 2 else 0) + (1 if len(roots2) == 2 else 0)
            n_total_fp = n_stable_states + n_crit_radii + (1 if (len(roots2) == 2 and not stable2) else 0)
            max_stable_droplets = 2 if stable2 else (1 if len(roots1) == 2 else 0)
        else:
            boundary_exists = False
            if R_single_bar is not None:
                _, dV2_b, _, _ = dV_pair(R_single_bar, 0.0)
                if dV2_b is not None and dV2_b < 0.0:
                    boundary_exists = True
            origin_stable = False
            dV1_0, dV2_0, _, _ = dV_pair(0.0, 0.0)
            if dV1_0 is not None and dV1_0 < 0.0 and dV2_0 < 0.0:
                origin_stable = True
            n_stable_states = (1 if stable2 else 0) + (2 if boundary_exists else 0) + (1 if origin_stable else 0)
            n_crit_radii = (1 if len(roots2) == 2 else 0) + (2 if len(roots1) == 2 else 0)
            n_total_fp = n_stable_states + n_crit_radii + (1 if (roots2 and not stable2) else 0)
            max_stable_droplets = 2 if stable2 else (1 if R_single_bar is not None else 0)

        return {
            "R_single_bar": R_single_bar,
            "R1_pair": R1_pair,
            "R2_pair": R2_pair,
            "V1_over_V2": V1_over_V2,
            "lam_pair": lam_pair,
            "n_stable_states": float(n_stable_states),
            "n_crit_radii": float(n_crit_radii),
            "n_total_fixed_points": float(n_total_fp),
            "max_stable_droplets": float(max_stable_droplets),
        }

    # Unequal cores Q_1 != Q_2: solve in log-radius coordinates (u1, u2) = (log(R1/R0), log(R2/R0))
    # so radii remain strictly positive and scale-independent across SI and dimensionless units
    V_est = max((V_s / psi) * max(phibar - k_b / k, 1e-9), 1e-30)
    frac1 = Q_1 / (Q_1 + Q_2)
    R1_guess = (3.0 * V_est * frac1 / (4.0 * math.pi)) ** (1.0 / 3.0)
    R2_guess = (3.0 * V_est * (1.0 - frac1) / (4.0 * math.pi)) ** (1.0 / 3.0)
    R_norm = (3.0 * V_est / (4.0 * math.pi)) ** (1.0 / 3.0)
    scale_rate = max(k_b * R_norm, 1e-30)

    def residual_u(u):
        r1 = R_norm * math.exp(max(min(float(u[0]), 20.0), -20.0))
        r2 = R_norm * math.exp(max(min(float(u[1]), 20.0), -20.0))
        dv1, dv2, _, _ = dV_pair(r1, r2)
        if dv1 is None or dv2 is None:
            return [-1e6, -1e6]
        return [
            dv1 / (4.0 * math.pi * r1**2 * scale_rate),
            dv2 / (4.0 * math.pi * r2**2 * scale_rate),
        ]

    sol = root(
        residual_u,
        [math.log(R1_guess / R_norm), math.log(R2_guess / R_norm)],
        tol=rtol,
    )
    if not sol.success:
        return {
            "R_single_bar": R_single_bar,
            "R1_pair": None,
            "R2_pair": None,
            "V1_over_V2": None,
            "lam_pair": None,
            "n_stable_states": 1.0 if R_single_bar is not None else 0.0,
            "n_crit_radii": 0.0,
            "n_total_fixed_points": 1.0 if R_single_bar is not None else 0.0,
            "max_stable_droplets": 1.0 if R_single_bar is not None else 0.0,
        }

    R1s = R_norm * math.exp(float(sol.x[0]))
    R2s = R_norm * math.exp(float(sol.x[1]))
    h = 1e-5 * min(R1s, R2s)
    J2 = np.zeros((2, 2), dtype=float)
    for j, (dr1, dr2) in enumerate([(h, 0.0), (0.0, h)]):
        dp1, dp2, _, _ = dV_pair(R1s + dr1, R2s + dr2)
        dm1, dm2, _, _ = dV_pair(R1s - dr1, R2s - dr2)
        J2[0, j] = (dp1 / (4.0 * math.pi * (R1s + dr1) ** 2) - dm1 / (4.0 * math.pi * (R1s - dr1) ** 2)) / (2.0 * h)
        J2[1, j] = (dp2 / (4.0 * math.pi * (R2s + dr2) ** 2) - dm2 / (4.0 * math.pi * (R2s - dr2) ** 2)) / (2.0 * h)
    evs = np.linalg.eigvals(J2)
    lam_max = float(np.max(evs.real))
    stable_pair = lam_max < 0.0

    n_boundary = 0
    if R_single_bar is not None:
        _, dV2_b, _, _ = dV_pair(R_single_bar, 0.0)
        if dV2_b is not None and dV2_b < 0.0:
            n_boundary += 1
    vals1_q2 = [dV_pair(0.0, float(r))[1] for r in grid]
    roots1_q2 = []
    for i in range(len(grid) - 1):
        if vals1_q2[i] is not None and vals1_q2[i + 1] is not None and vals1_q2[i] * vals1_q2[i + 1] < 0.0:
            roots1_q2.append(
                float(brentq(lambda r: dV_pair(0.0, r)[1], float(grid[i]), float(grid[i + 1]), xtol=rtol * R_lo, rtol=rtol))
            )
    if roots1_q2:
        dV1_b, _, _, _ = dV_pair(0.0, roots1_q2[-1])
        if dV1_b is not None and dV1_b < 0.0:
            n_boundary += 1

    n_total_fp = (1 if stable_pair else 0) + n_boundary
    return {
        "R_single_bar": R_single_bar,
        "R1_pair": R1s,
        "R2_pair": R2s,
        "V1_over_V2": (R1s / R2s) ** 3,
        "lam_pair": lam_max,
        "n_stable_states": float(n_total_fp),
        "n_crit_radii": 0.0,
        "n_total_fixed_points": float(n_total_fp),
        "max_stable_droplets": 2.0 if stable_pair else 1.0,
    }
