def run(inputs, settings):
    return {"c_onset": inputs["c_KCl"] ** inputs["N"]}


def symbolic():
    import sympy as sp

    c_KCl, N = sp.symbols("c_KCl N", positive=True)
    return {"c_onset": c_KCl ** N}
