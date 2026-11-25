def run(inputs, settings):
    # (uN/m) / um = 1 N/m^2
    return {"DeltaP_column": inputs["sigma"] / inputs["R_f"]}


def symbolic():
    import sympy as sp
    sigma, R_f = sp.symbols("sigma R_f", positive=True)
    return {"DeltaP_column": sigma / R_f}
