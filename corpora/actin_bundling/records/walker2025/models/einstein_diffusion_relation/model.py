def run(inputs, settings):
    # mu in um/(pN*s), kBT in pN*nm (1 pN*nm = 1e-3 pN*um) -> D in um^2/s
    return {"D": float(inputs["mu"]) * float(inputs["kBT"]) * 1e-3}


def symbolic():
    import sympy as sp

    mu, kBT = sp.symbols("mu kBT", positive=True)
    return {"D": mu * kBT * sp.Rational(1, 1000)}
