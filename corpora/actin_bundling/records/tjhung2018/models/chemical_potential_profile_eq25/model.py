def run(inputs, settings):
    A = float(inputs["A"])
    mu_I = float(inputs["mu_I"])
    phi_s = float(inputs["phi_s"])
    R = float(inputs["R"])
    r = float(inputs["r"])
    mu_inf = A * (phi_s**3 - phi_s)
    return {"mu": R * (mu_I - mu_inf) / r + mu_inf}


def symbolic():
    import sympy as sp

    A, mu_I, phi_s, R, r = sp.symbols("A mu_I phi_s R r", real=True)
    mu_inf = A * (phi_s**3 - phi_s)
    return {"mu": R * (mu_I - mu_inf) / r + mu_inf}
