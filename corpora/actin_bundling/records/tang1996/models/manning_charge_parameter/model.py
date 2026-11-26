def run(inputs, settings):
    return {"zeta": inputs["lambda_B"] / inputs["b"]}


def symbolic():
    import sympy as sp

    lambda_B, b = sp.symbols("lambda_B b", positive=True)
    return {"zeta": lambda_B / b}
