def run(inputs, settings):
    return {"dE": inputs["dE_kBT"] * inputs["kBT"]}


def symbolic():
    import sympy as sp
    dE_kBT, kBT = sp.symbols("dE_kBT kBT", positive=True)
    return {"dE": dE_kBT * kBT}
