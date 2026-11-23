import math


def run(inputs, settings):
    R = inputs["R"]
    gamma = inputs["gamma"]
    mu_perp_bar = inputs["mu_perp_bar"]
    K3_bar = inputs["K3_bar"]
    epsilon = inputs["epsilon"]
    Sigma = inputs["Sigma"]
    lam3_sq = K3_bar / mu_perp_bar
    om0 = (
        (gamma / (4.0 * mu_perp_bar)) ** (1.0 / 3.0)
        * (R ** (-2.0 / 3.0))
        * ((R**2 + lam3_sq) ** (-1.0 / 3.0))
    )
    e_elast = mu_perp_bar * (om0 * R) ** 4 + K3_bar * (om0**4) * (R**2) - gamma * om0
    elast = math.pi * (R**2) * e_elast
    total = elast - math.pi * (R**2) * epsilon + 2.0 * math.pi * R * Sigma
    return {"dF_over_L": total, "dF_elast_over_L": elast}


def symbolic():
    import sympy as sp
    R, gamma, mu_perp_bar, K3_bar, epsilon, Sigma = sp.symbols(
        "R gamma mu_perp_bar K3_bar epsilon Sigma", positive=True
    )
    lam3_sq = K3_bar / mu_perp_bar
    om0 = (
        (gamma / (4 * mu_perp_bar)) ** sp.Rational(1, 3)
        * R ** (-sp.Rational(2, 3))
        * (R**2 + lam3_sq) ** (-sp.Rational(1, 3))
    )
    e_elast = mu_perp_bar * (om0 * R) ** 4 + K3_bar * om0**4 * R**2 - gamma * om0
    elast = sp.simplify(sp.pi * R**2 * e_elast)
    total = elast - sp.pi * R**2 * epsilon + 2 * sp.pi * R * Sigma
    return {"dF_over_L": total, "dF_elast_over_L": elast}
