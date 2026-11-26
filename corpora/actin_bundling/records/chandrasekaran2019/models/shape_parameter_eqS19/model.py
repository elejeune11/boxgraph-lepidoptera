def run(inputs, settings):
    l1, l2, l3 = inputs["lambda_1"], inputs["lambda_2"], inputs["lambda_3"]
    return {"Sh": l1 / (l1 + l2 + l3)}


def symbolic():
    import sympy as sp
    lambda_1, lambda_2, lambda_3 = sp.symbols("lambda_1 lambda_2 lambda_3", positive=True)
    return {"Sh": lambda_1 / (lambda_1 + lambda_2 + lambda_3)}
