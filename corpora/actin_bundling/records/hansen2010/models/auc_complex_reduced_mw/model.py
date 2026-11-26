def run(inputs, settings):
    return {"sigma_i": inputs["sigma_o"] + inputs["i"] * inputs["sigma_a"]}


def symbolic():
    import sympy as sp

    sigma_o, sigma_a, i = sp.symbols("sigma_o sigma_a i", positive=True)
    return {"sigma_i": sigma_o + i * sigma_a}
