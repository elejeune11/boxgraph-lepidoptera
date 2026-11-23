def run(inputs, settings):
    k2 = inputs["k"] ** 2
    lam = -k2 * inputs["D"] + (
        k2 * inputs["c0"] * inputs["zdm0"] * inputs["df_dc"]
        / (inputs["gamma"] + inputs["eta"] * k2)
    )
    return {"lam": lam}


def symbolic():
    import sympy as sp

    k, D, c0, zdm0, df_dc, gamma, eta = sp.symbols("k D c0 zdm0 df_dc gamma eta", positive=True)
    return {"lam": -k**2 * D + k**2 * c0 * zdm0 * df_dc / (gamma + eta * k**2)}
