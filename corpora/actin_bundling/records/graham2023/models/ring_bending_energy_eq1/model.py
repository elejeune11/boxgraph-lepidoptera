def run(inputs, settings):
    # k_b in pN*nm^2; L, R, u_m in um -> 1/um = 1e-3 / nm -> E_bend in pN*nm
    k_b = inputs["k_b"]
    L = inputs["L"]
    R = inputs["R"]
    m = inputs["m"]
    u_m = inputs["u_m"]
    return {"E_bend": 0.5 * k_b * 1e-3 * (L / (R ** 2) + (L / (R ** 4)) * (2.0 + m ** 4) * (u_m ** 2))}


def symbolic():
    import sympy as sp

    k_b, L, R, m, u_m = sp.symbols("k_b L R m u_m", positive=True)
    return {"E_bend": sp.Rational(1, 2000) * k_b * (L / R ** 2 + (L / R ** 4) * (2 + m ** 4) * u_m ** 2)}
