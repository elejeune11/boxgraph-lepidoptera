def run(inputs, settings):
    k_b = inputs["k_b"]
    L = inputs["L"]
    R = inputs["R"]
    m = inputs["m"]
    u_m = inputs["u_m"]
    coeff = 0.5 * (1.0 + 6.0 * (m ** 2) + (m ** 4))
    return {"E_bend": 0.5 * k_b * 1e-3 * (L / (R ** 2) + (L / (R ** 4)) * coeff * (u_m ** 2))}


def symbolic():
    import sympy as sp

    k_b, L, R, m, u_m = sp.symbols("k_b L R m u_m", positive=True)
    coeff = sp.Rational(1, 2) * (1 + 6 * m ** 2 + m ** 4)
    return {"E_bend": sp.Rational(1, 2000) * k_b * (L / R ** 2 + (L / R ** 4) * coeff * u_m ** 2)}
