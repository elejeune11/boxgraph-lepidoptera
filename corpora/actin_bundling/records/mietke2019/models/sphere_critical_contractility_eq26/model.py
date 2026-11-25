def run(inputs, settings):
    tau = inputs["eta_b"] / inputs["gamma"]
    return {"alpha_s_star": tau * (2.0 * inputs["D"] / (inputs["R0"] ** 2) + inputs["k"])}


def symbolic():
    import sympy as sp

    eta_b, gamma, D, R0, k = sp.symbols("eta_b gamma D R0 k", nonnegative=True)
    return {"alpha_s_star": (eta_b / gamma) * (2 * D / R0**2 + k)}
