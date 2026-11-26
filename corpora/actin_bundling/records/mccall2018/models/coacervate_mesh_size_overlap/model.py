def run(inputs, settings):
    a = inputs["a"]
    N = inputs["N"]
    Phi_pp = inputs["Phi_pp"]
    R_F = a * (N ** (3.0 / 5.0))
    Phi_star = (N / (R_F ** 3)) * (a ** 3)
    return {"xi": R_F * ((Phi_star / Phi_pp) ** (3.0 / 4.0))}


def symbolic():
    import sympy as sp

    a, N, Phi_pp = sp.symbols("a N Phi_pp", positive=True)
    R_F = a * N ** sp.Rational(3, 5)
    Phi_star = (N / R_F ** 3) * a ** 3
    return {"xi": R_F * (Phi_star / Phi_pp) ** sp.Rational(3, 4)}
