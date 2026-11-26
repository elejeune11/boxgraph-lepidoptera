def run(inputs, settings):
    return {"W_bar": inputs["X_1"] * inputs["W"]}


def symbolic():
    import sympy as sp
    X_1, W = sp.symbols("X_1 W", real=True)
    return {"W_bar": X_1 * W}
