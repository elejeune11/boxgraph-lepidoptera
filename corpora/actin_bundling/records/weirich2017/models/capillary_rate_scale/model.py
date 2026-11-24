def run(inputs, settings):
    # gamma_I [1e-9 N/m] / (eta [N*s/m^2] * R_0 [1e-9 m]) = gamma_I / (eta * R_0) [1/s]
    return {"rate_0": inputs["gamma_I"] / (inputs["eta"] * inputs["R_0"])}


def symbolic():
    import sympy as sp

    gamma_I, eta, R_0 = sp.symbols("gamma_I eta R_0", positive=True)
    return {"rate_0": gamma_I / (eta * R_0)}
