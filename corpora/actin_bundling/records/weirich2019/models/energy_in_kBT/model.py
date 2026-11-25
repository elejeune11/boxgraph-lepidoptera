def run(inputs, settings):
    return {"F0_kBT": inputs["F_0"] / inputs["kBT"]}


def symbolic():
    import sympy as sp
    F_0, kBT = sp.symbols("F_0 kBT", positive=True)
    return {"F0_kBT": F_0 / kBT}
