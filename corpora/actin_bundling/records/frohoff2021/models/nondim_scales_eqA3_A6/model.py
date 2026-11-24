import math


def run(inputs, settings):
    L = inputs["L"]
    L0 = inputs["L0"]
    Q1 = inputs["Q1"]
    Q2 = inputs["Q2"]
    kappa1 = inputs["kappa1"]
    kappa2 = inputs["kappa2"]
    zeta1 = inputs["zeta1"]
    zeta2 = inputs["zeta2"]
    a1 = inputs["a1"]
    a2 = inputs["a2"]
    b1 = inputs["b1"]
    b2 = inputs["b2"]
    gamma1 = inputs["gamma1"]
    gamma2 = inputs["gamma2"]

    tau = (L0**4) / (Q1 * kappa1)
    ell = L / L0
    Q = Q2 / Q1
    kappa = kappa2 / kappa1
    phi1_hat = math.sqrt(kappa1 / (zeta1 * b1 * L0**2))
    phi2_hat = math.sqrt(kappa1 / (zeta2 * b2 * L0**2))
    a_tilde_1 = a1 / (b1 * phi1_hat**2)
    a_tilde_2 = a2 / (b2 * phi2_hat**2)
    a = a_tilde_1
    a_Delta = a_tilde_2 - a_tilde_1
    gamma_tilde_1 = (phi2_hat * gamma1 * L0**2) / (phi1_hat * kappa1)
    gamma_tilde_2 = (phi1_hat * gamma2 * L0**2) / (phi2_hat * kappa1)
    rho = 0.5 * (gamma_tilde_1 + gamma_tilde_2)
    alpha = 0.5 * (gamma_tilde_1 - gamma_tilde_2)
    return {
        "tau": tau,
        "ell": ell,
        "Q": Q,
        "kappa": kappa,
        "phi1_hat": phi1_hat,
        "phi2_hat": phi2_hat,
        "a": a,
        "a_Delta": a_Delta,
        "rho": rho,
        "alpha": alpha,
    }


def symbolic():
    import sympy as sp

    L, L0, Q1, Q2, kappa1, kappa2 = sp.symbols("L L0 Q1 Q2 kappa1 kappa2", positive=True)
    zeta1, zeta2, b1, b2 = sp.symbols("zeta1 zeta2 b1 b2", positive=True)
    a1, a2, gamma1, gamma2 = sp.symbols("a1 a2 gamma1 gamma2", real=True)
    phi1_hat = sp.sqrt(kappa1 / (zeta1 * b1 * L0**2))
    phi2_hat = sp.sqrt(kappa1 / (zeta2 * b2 * L0**2))
    a_tilde_1 = a1 / (b1 * phi1_hat**2)
    a_tilde_2 = a2 / (b2 * phi2_hat**2)
    gt1 = (phi2_hat * gamma1 * L0**2) / (phi1_hat * kappa1)
    gt2 = (phi1_hat * gamma2 * L0**2) / (phi2_hat * kappa1)
    return {
        "tau": L0**4 / (Q1 * kappa1),
        "ell": L / L0,
        "Q": Q2 / Q1,
        "kappa": kappa2 / kappa1,
        "phi1_hat": phi1_hat,
        "phi2_hat": phi2_hat,
        "a": a_tilde_1,
        "a_Delta": a_tilde_2 - a_tilde_1,
        "rho": (gt1 + gt2) / 2,
        "alpha": (gt1 - gt2) / 2,
    }
