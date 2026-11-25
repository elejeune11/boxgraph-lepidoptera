def run(inputs, settings):
    return {"alpha_c_star": inputs["D_tilde"] + inputs["k_tilde"]}


def symbolic():
    import sympy as sp

    D_tilde, k_tilde = sp.symbols("D_tilde k_tilde", nonnegative=True)
    return {"alpha_c_star": D_tilde + k_tilde}
