def run(inputs, settings):
    return {"rho_v": inputs["rho"] / inputs["L_p"]}


def symbolic():
    import sympy as sp

    rho, L_p = sp.symbols("rho L_p", positive=True)
    return {"rho_v": rho / L_p}
