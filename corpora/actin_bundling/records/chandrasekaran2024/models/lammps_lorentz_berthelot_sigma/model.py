def run(inputs, settings):
    return {"sigma_12": 0.5 * (inputs["sigma_1"] + inputs["sigma_2"])}


def symbolic():
    import sympy as sp

    sigma_1, sigma_2 = sp.symbols("sigma_1 sigma_2", positive=True)
    return {"sigma_12": (sigma_1 + sigma_2) / 2}
