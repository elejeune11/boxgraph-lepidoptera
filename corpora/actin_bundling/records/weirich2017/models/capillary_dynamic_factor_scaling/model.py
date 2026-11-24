def run(inputs, settings):
    # gamma_I [1e-9 N/m] / (eta [N*s/m^2] * R_0 [1e-9 m]) = [1/s]
    return {"Lambda": (inputs["gamma_I"] / (inputs["eta"] * inputs["R_0"])) * (inputs["k"] ** 2)}


def symbolic():
    import sympy as sp

    gamma_I, eta, R_0, k = sp.symbols("gamma_I eta R_0 k", positive=True)
    return {"Lambda": (gamma_I / (eta * R_0)) * k**2}
