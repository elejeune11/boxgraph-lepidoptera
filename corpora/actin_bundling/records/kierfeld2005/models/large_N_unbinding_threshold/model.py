def run(inputs, settings):
    return {"W_bar_c": 2.0 * inputs["W_bar_c_pair"] / inputs["q"]}


def symbolic():
    import sympy as sp
    W_bar_c_pair, q = sp.symbols("W_bar_c_pair q", positive=True)
    return {"W_bar_c": 2 * W_bar_c_pair / q}
