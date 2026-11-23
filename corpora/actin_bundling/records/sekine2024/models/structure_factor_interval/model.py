"""Structure interval = 2*pi/mu = 1/(mu/(2*pi)) converted from um to nm (Methods p. 11 & Table 4)."""


def run(inputs, settings):
    # 1 / (1/um) = um; multiply by 1000 for nm
    return {"interval": 1000.0 / inputs["mu_over_2pi"]}


def symbolic():
    import sympy as sp

    mu_over_2pi = sp.Symbol("mu_over_2pi", positive=True)
    return {"interval": sp.Integer(1000) / mu_over_2pi}
