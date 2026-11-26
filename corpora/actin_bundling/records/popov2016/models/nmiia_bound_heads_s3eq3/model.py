def run(inputs, settings):
    return {"N_bound": inputs["rho"] * inputs["N_total"] + inputs["gamma"] * inputs["F_ext"]}

def symbolic():
    import sympy as sp
    rho, N_total, gamma, F_ext = sp.symbols("rho N_total gamma F_ext", nonnegative=True)
    return {"N_bound": rho * N_total + gamma * F_ext}
