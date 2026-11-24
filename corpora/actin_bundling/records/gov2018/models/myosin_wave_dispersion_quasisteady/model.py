import math


def run(inputs, settings):
    q = inputs["q"]
    kappa = inputs["kappa"]
    sigma = inputs["sigma"]
    kBT = inputs["kBT"]
    Lambda = inputs["Lambda"]
    H0 = inputs["H0"]
    ns = inputs["ns"]
    n0 = inputs["n0"]
    A = inputs["A"]
    A_star = inputs["A_star"]
    k_on = inputs["k_on"]
    k_off = inputs["k_off"]

    q2 = q * q
    denom = ns * (kappa * q2 + sigma)
    # Exact xi -> 0 quadratic: omega^2 + b*omega + c0 = 0 (convex H0 < 0 -> -H0 = |H0| > 0)
    b = (
        A * H0 * Lambda * kappa * n0 * q2
        + (H0 ** 2) * Lambda * kappa * n0 * q2 * sigma
        + Lambda * kBT * kappa * ns * q2 * q2
        + Lambda * kBT * ns * q2 * sigma
        + k_off * denom
    ) / denom
    c0 = (
        Lambda
        * q2
        * (
            A * H0 * kappa * k_off * n0
            - A_star * H0 * kappa * k_on * n0
            + (H0 ** 2) * kappa * k_off * n0 * sigma
            + kBT * kappa * k_off * ns * q2
            + kBT * k_off * ns * sigma
        )
        / denom
    )
    disc = b * b - 4.0 * c0
    if disc >= 0.0:
        omega_real = 0.5 * (-b + math.sqrt(disc))
        omega_imag = 0.0
    else:
        omega_real = -0.5 * b
        omega_imag = 0.5 * math.sqrt(-disc)

    H0_mag = abs(H0)
    c_ana = math.sqrt((kappa * Lambda * H0_mag * n0 * k_on * A_star) / (sigma * ns)) if A_star > 0.0 else 0.0
    q_omega_approx = k_off / (2.0 * c_ana) if c_ana > 0.0 else None
    return {
        "omega_real": omega_real,
        "omega_imag": omega_imag,
        "c_ana": c_ana,
        "q_omega_approx": q_omega_approx,
    }


def symbolic():
    import sympy as sp

    q, kappa, sigma, kBT, Lambda, H0, ns, n0, A, A_star, k_on, k_off = sp.symbols(
        "q kappa sigma kBT Lambda H0 ns n0 A A_star k_on k_off", positive=True
    )
    # Note: in positive-symbol form, H0 represents |H0| (for convex CMC H0_signed = -H0)
    c_ana = sp.sqrt((kappa * Lambda * H0 * n0 * k_on * A_star) / (sigma * ns))
    return {
        "c_ana": c_ana,
        "q_omega_approx": k_off / (2 * c_ana),
    }
