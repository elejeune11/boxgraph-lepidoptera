def run(inputs, settings):
    return {"R_0": 1.0 / inputs["rho"]}


def symbolic():
    import sympy as sp

    rho = sp.symbols("rho", positive=True)
    return {"R_0": 1 / rho}
