def run(inputs, settings):
    return {"X_1c": inputs["W_bar_c"] / inputs["W"]}


def symbolic():
    import sympy as sp
    W_bar_c, W = sp.symbols("W_bar_c W", real=True)
    return {"X_1c": W_bar_c / W}
