def run(inputs, settings):
    return {"L": inputs["d_mon"] * (inputs["c_actin"] * 1000.0) / inputs["c_cp"]}


def symbolic():
    import sympy as sp

    d_mon, c_actin, c_cp = sp.symbols("d_mon c_actin c_cp", positive=True)
    return {"L": 1000 * d_mon * c_actin / c_cp}
