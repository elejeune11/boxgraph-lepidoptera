def run(inputs, settings):
    return {"f_neut": inputs["theta_18"] * inputs["N"]}


def symbolic():
    import sympy as sp

    theta_18, N = sp.symbols("theta_18 N", positive=True)
    return {"f_neut": theta_18 * N}
