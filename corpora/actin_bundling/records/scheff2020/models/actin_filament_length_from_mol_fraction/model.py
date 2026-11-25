def run(inputs, settings):
    # f_cp is in percent (1% = 0.01)
    return {"L": 100.0 * inputs["d_mon"] / inputs["f_cp"]}


def symbolic():
    import sympy as sp

    d_mon, f_cp = sp.symbols("d_mon f_cp", positive=True)
    return {"L": 100 * d_mon / f_cp}
