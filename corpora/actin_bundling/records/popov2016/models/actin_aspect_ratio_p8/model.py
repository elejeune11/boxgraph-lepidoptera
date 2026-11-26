def run(inputs, settings):
    return {"aspect_ratio": inputs["l_p"] / inputs["sigma_0"]}

def symbolic():
    import sympy as sp
    l_p, sigma_0 = sp.symbols("l_p sigma_0", positive=True)
    return {"aspect_ratio": l_p / sigma_0}
