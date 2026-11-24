def run(inputs, settings):
    return {"k_bend_thermal": inputs["k_bend_phys"] / inputs["kBT"]}


def symbolic():
    import sympy as sp

    k_bend_phys, kBT = sp.symbols("k_bend_phys kBT", positive=True)
    return {"k_bend_thermal": k_bend_phys / kBT}
