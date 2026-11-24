def run(inputs, settings):
    xi = inputs["xi"]
    return {"c": (0.3 / xi) ** 2}


def symbolic():
    import sympy as sp

    xi = sp.Symbol("xi", positive=True)
    return {"c": (sp.Rational(3, 10) / xi) ** 2}
