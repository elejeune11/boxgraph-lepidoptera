def run(inputs, settings):
    # 1 uM = 1000 nM
    return {"N_m": (inputs["c_actin"] * 1000.0) / inputs["c_cp"]}


def symbolic():
    import sympy as sp

    c_actin, c_cp = sp.symbols("c_actin c_cp", positive=True)
    return {"N_m": 1000 * c_actin / c_cp}
