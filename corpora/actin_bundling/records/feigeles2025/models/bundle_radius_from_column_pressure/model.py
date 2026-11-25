def run(inputs, settings):
    # (uN/m) / (N/m^2) = 1e-6 m = 1 um
    return {"R_f": inputs["sigma"] / inputs["DeltaP_column"]}


def symbolic():
    import sympy as sp
    sigma, DeltaP_column = sp.symbols("sigma DeltaP_column", positive=True)
    return {"R_f": sigma / DeltaP_column}
