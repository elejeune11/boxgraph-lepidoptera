def run(inputs, settings):
    # 1 (pN/um) / um = 1 pN/um^2 = 1e-12 N / 1e-12 m^2 = 1 Pa
    return {"p0_s": 2.0 * (inputs["gamma"] + inputs["xi"] * inputs["f_c0"]) / inputs["R0"]}


def symbolic():
    import sympy as sp

    gamma, xi, f_c0, R0 = sp.symbols("gamma xi f_c0 R0", positive=True)
    return {"p0_s": 2 * (gamma + xi * f_c0) / R0}
