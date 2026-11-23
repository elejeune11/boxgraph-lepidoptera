def run(inputs, settings):
    return {"pct_increase": 100.0 * (inputs["lambda_mutant"] - inputs["lambda_wt"]) / inputs["lambda_wt"]}


def symbolic():
    import sympy as sp

    lambda_mutant, lambda_wt = sp.symbols("lambda_mutant lambda_wt", positive=True)
    return {"pct_increase": 100 * (lambda_mutant - lambda_wt) / lambda_wt}
