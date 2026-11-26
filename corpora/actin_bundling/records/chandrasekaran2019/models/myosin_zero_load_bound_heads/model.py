def run(inputs, settings):
    return {"N_0": inputs["rho"] * inputs["N_total"]}


def symbolic():
    import sympy as sp
    rho, N_total = sp.symbols("rho N_total", positive=True)
    return {"N_0": rho * N_total}
