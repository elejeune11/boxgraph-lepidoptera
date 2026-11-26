"""Exact Fourier mode-sum three-point bending stiffness from Appendix Eqs. (A6)-(A9)."""

import math


def run(inputs, settings):
    kappa_f = inputs["kappa_f"]
    N = inputs["N"]
    alpha = inputs["alpha"]
    chi2 = inputs["chi2"]
    bc = settings.get("bc", "pinned")
    n_modes = int(settings.get("n_modes", 200))

    # For midpoint point load F_{L/2}:
    # - pinned (sin(j*pi*x/L), odd j = 1, 3, 5, ...): q_j*L = j*pi, sum_{odd j} 2 / (j*pi)^4 = 1/48 (beta = 48)
    # - clamped (cos(2*j*pi*x/L) - 1, odd j = 1, 3, 5, ...): q_j*L = 2*j*pi, sum_{odd j} 8 / (2*j*pi)^4 = 1/192 (beta = 192)
    q_mult = math.pi if bc == "pinned" else 2.0 * math.pi
    beta = 48.0 if bc == "pinned" else 192.0
    pref = 2.0 if bc == "pinned" else 8.0

    # Dimensionless compliance: |r_perp(L/2)| * kappa_f / (F_{L/2} * L^3) = sum_{k=0}^{n_modes-1} pref / ((q_j*L)^4 * kappa_B*(q_j))
    sqrt_N = math.sqrt(N)
    comp = 0.0
    for k in range(n_modes):
        j = 2 * k + 1
        qL = j * q_mult
        c_q = (qL * qL) / 12.0
        kb_star_q = N * (1.0 + chi2 * (N - 1.0) / (1.0 + c_q * (N + sqrt_N) / alpha))
        comp += pref / ((qL ** 4) * kb_star_q)

    kappa_B_star = 1.0 / (beta * comp)
    return {"kappa_B": kappa_f * kappa_B_star, "kappa_B_star": kappa_B_star}
