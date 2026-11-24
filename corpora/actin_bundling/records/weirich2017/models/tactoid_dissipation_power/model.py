def run(inputs, settings):
    # V [um^3 = 1e-18 m^3] * eta [Pa*s = J*s/m^3] * (1/s)^2 = 1e-18 J/s = 1e-18 W
    V_m3 = inputs["V"] * 1e-18
    term_R = inputs["eta_R"] * (inputs["dR_dt"] / inputs["R"]) ** 2
    term_a = inputs["eta_alpha"] * (inputs["dalpha_dt"] / inputs["alpha"]) ** 2
    return {"D": 0.5 * V_m3 * (term_R + term_a)}


def symbolic():
    import sympy as sp

    V, eta_R, eta_alpha, R, dR_dt, alpha, dalpha_dt = sp.symbols(
        "V eta_R eta_alpha R dR_dt alpha dalpha_dt", positive=True
    )
    return {
        "D": sp.Rational(1, 2)
        * sp.Rational(1, 10**18)
        * V
        * (eta_R * (dR_dt / R) ** 2 + eta_alpha * (dalpha_dt / alpha) ** 2)
    }
