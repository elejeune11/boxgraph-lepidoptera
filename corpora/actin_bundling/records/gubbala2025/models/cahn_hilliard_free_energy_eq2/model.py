def run(inputs, settings):
    phi = float(inputs["phi"])
    lam = float(inputs["lam"])
    g = float(inputs["grad_phi"])
    f_phi = 0.25 * phi**4 - 0.5 * phi**2 + 0.5 * lam * g**2
    df_dphi = phi**3 - phi
    return {"f_phi": f_phi, "df_dphi": df_dphi}


def symbolic():
    import sympy as sp

    phi, lam, grad_phi = sp.symbols("phi lam grad_phi", real=True)
    return {
        "f_phi": phi**4 / 4 - phi**2 / 2 + lam * grad_phi**2 / 2,
        "df_dphi": phi**3 - phi,
    }
