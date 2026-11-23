def run(inputs, settings):
    return {"lambda_R": inputs["lambda_R_phys"] / inputs["R_g"]}


def symbolic():
    import sympy as sp
    lambda_R_phys, R_g = sp.symbols("lambda_R_phys R_g", positive=True)
    return {"lambda_R": lambda_R_phys / R_g}
