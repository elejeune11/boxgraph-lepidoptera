def run(inputs, settings):
    return {"D": inputs["R_star"] ** 0.3}


def symbolic():
    import sympy as sp

    R_star = sp.Symbol("R_star", positive=True)
    return {"D": R_star ** sp.Rational(3, 10)}
