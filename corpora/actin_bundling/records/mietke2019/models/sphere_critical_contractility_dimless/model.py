def run(inputs, settings):
    return {"alpha_s_star": 2.0 * inputs["D_tilde"] + inputs["k_tilde"]}


def symbolic():
    import sympy as sp

    D_tilde, k_tilde = sp.symbols("D_tilde k_tilde", nonnegative=True)
    return {"alpha_s_star": 2 * D_tilde + k_tilde}
