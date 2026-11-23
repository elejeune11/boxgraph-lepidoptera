def run(inputs, settings):
    return {"fold": inputs["lambda_mutant"] / inputs["lambda_wt"]}


def symbolic():
    import sympy as sp

    lambda_mutant, lambda_wt = sp.symbols("lambda_mutant lambda_wt", positive=True)
    return {"fold": lambda_mutant / lambda_wt}
