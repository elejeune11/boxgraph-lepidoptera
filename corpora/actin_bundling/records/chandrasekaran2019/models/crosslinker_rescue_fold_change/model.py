def run(inputs, settings):
    return {"fold": inputs["alpha_high"] / inputs["alpha_low"]}


def symbolic():
    import sympy as sp
    alpha_high, alpha_low = sp.symbols("alpha_high alpha_low", positive=True)
    return {"fold": alpha_high / alpha_low}
