def run(inputs, settings):
    return {"N_bound": inputs["N_0"] + inputs["gamma"] * inputs["F_ext"]}


def symbolic():
    import sympy as sp
    N_0, gamma, F_ext = sp.symbols("N_0 gamma F_ext", nonnegative=True)
    return {"N_bound": N_0 + gamma * F_ext}
