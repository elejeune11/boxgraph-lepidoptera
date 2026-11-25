def run(inputs, settings):
    w_pN_um = inputs["w"] * 1e-3          # nN/m -> pN/um
    gamma_pN_um = inputs["gamma"] * 1e-3  # nN/m -> pN/um
    K = inputs["K"]                       # pN
    a = inputs["a"]                       # um
    W_sph = w_pN_um / gamma_pN_um
    W_el = w_pN_um * (a ** 2) * gamma_pN_um / (K ** 2)
    return {"W_sph": W_sph, "W_el": W_el}


def symbolic():
    import sympy as sp
    w, gamma, K, a = sp.symbols("w gamma K a", positive=True)
    w_pN = w * sp.Rational(1, 1000)
    gamma_pN = gamma * sp.Rational(1, 1000)
    return {
        "W_sph": w / gamma,
        "W_el": w_pN * a ** 2 * gamma_pN / K ** 2,
    }
