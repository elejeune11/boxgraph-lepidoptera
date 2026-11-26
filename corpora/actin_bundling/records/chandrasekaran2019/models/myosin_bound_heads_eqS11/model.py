def run(inputs, settings):
    # As printed on p. 34, gamma = 0.05 lacks its 1/pN unit; evaluated numerically with F_ext in pN.
    return {"N_bound": inputs["N_0"] + inputs["gamma"] * inputs["F_ext"]}


def symbolic():
    import sympy as sp
    N_0, gamma, F_ext = sp.symbols("N_0 gamma F_ext", nonnegative=True)
    return {"N_bound": N_0 + gamma * F_ext}
