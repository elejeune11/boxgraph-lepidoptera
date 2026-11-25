def run(inputs, settings):
    tau = inputs["eta_b"] / inputs["gamma"]
    return {"alpha_c_star": tau * (inputs["D"] / (inputs["r0"] ** 2) + inputs["k"])}


def symbolic():
    import sympy as sp

    eta_b, gamma, D, r0, k = sp.symbols("eta_b gamma D r0 k", nonnegative=True)
    return {"alpha_c_star": (eta_b / gamma) * (D / r0**2 + k)}
