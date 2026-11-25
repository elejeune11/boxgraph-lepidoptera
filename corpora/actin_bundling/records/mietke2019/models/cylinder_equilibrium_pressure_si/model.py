def run(inputs, settings):
    # kappa [pN*um] / r0^2 [um^2] = pN/um; dividing by r0 [um] gives pN/um^2 = Pa
    r0 = inputs["r0"]
    return {"p0_c": (inputs["gamma"] + inputs["xi"] * inputs["f_c0"] - inputs["kappa"] / (r0 ** 2)) / r0}


def symbolic():
    import sympy as sp

    gamma, xi, f_c0, kappa, r0 = sp.symbols("gamma xi f_c0 kappa r0", nonnegative=True)
    return {"p0_c": (gamma + xi * f_c0 - kappa / r0**2) / r0}
