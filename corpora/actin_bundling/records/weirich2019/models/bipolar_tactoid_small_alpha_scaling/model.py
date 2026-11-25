def run(inputs, settings):
    K = inputs["K"]
    gamma_pN = inputs["gamma"] * 1e-3  # nN/m -> pN/um
    R = inputs["R"]
    alpha = inputs["alpha"]
    return {
        "V": (R ** 3) * (alpha ** 5),
        "F_bare": ((K ** 2) / gamma_pN) * (alpha ** 3),
    }


def symbolic():
    import sympy as sp
    K, gamma, R, alpha = sp.symbols("K gamma R alpha", positive=True)
    gamma_pN = gamma * sp.Rational(1, 1000)
    return {
        "V": R ** 3 * alpha ** 5,
        "F_bare": (K ** 2 / gamma_pN) * alpha ** 3,
    }
