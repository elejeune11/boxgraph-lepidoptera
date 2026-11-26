def run(inputs, settings):
    return {"ratio": inputs["mu"] / inputs["mu_water"]}


def symbolic():
    import sympy as sp

    mu, mu_water = sp.symbols("mu mu_water", positive=True)
    return {"ratio": mu / mu_water}
