def run(inputs, settings):
    A_actin = inputs["A_total"] - inputs["A_Lo"]
    rho = inputs["I_actin"] / A_actin
    return {"A_actin": A_actin, "rho": rho}


def symbolic():
    import sympy as sp

    A_total, A_Lo, I_actin = sp.symbols("A_total A_Lo I_actin", positive=True)
    A_actin = A_total - A_Lo
    return {"A_actin": A_actin, "rho": I_actin / A_actin}
