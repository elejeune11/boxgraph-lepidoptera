def run(inputs, settings):
    return {"ratio": inputs["c_dense"] / inputs["c_dilute"]}


def symbolic():
    import sympy as sp

    c_dense, c_dilute = sp.symbols("c_dense c_dilute", positive=True)
    return {"ratio": c_dense / c_dilute}
