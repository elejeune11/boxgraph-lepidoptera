def run(inputs, settings):
    q2 = inputs["k_l"] ** 2
    c0 = inputs["c0"]
    alpha = inputs["Pe"] * c0 / ((1.0 + c0) ** 2)
    return {"lam_dimless": -q2 * (1.0 - alpha / (1.0 + q2))}


def symbolic():
    import sympy as sp

    k_l, Pe, c0 = sp.symbols("k_l Pe c0", positive=True)
    return {"lam_dimless": -(k_l**2) * (1 - Pe * c0 / ((1 + c0) ** 2 * (1 + k_l**2)))}
