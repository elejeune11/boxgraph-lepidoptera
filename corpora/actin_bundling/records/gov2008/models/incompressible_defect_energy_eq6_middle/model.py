def run(inputs, settings):
    K = inputs["K"]
    rho_v = inputs["rho_v"]
    L_p = inputs["L_p"]
    a = inputs["a"]
    R = inputs["R"]
    return {"E_slice": K * rho_v * L_p * ((a / (L_p * L_p)) ** 2) * (R * R) * L_p}


def symbolic():
    import sympy as sp

    K, rho_v, L_p, a, R = sp.symbols("K rho_v L_p a R", positive=True)
    return {"E_slice": K * rho_v * L_p * (a / L_p**2) ** 2 * R**2 * L_p}
