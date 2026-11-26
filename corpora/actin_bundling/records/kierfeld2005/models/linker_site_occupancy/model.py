def run(inputs, settings):
    kc = inputs["K"] * inputs["c_x"]
    return {"X_1": kc / (1.0 + kc)}


def symbolic():
    import sympy as sp
    K, c_x = sp.symbols("K c_x", positive=True)
    return {"X_1": K * c_x / (1 + K * c_x)}
