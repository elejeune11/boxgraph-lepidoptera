def run(inputs, settings):
    return {"mu_par_bar": 0.25 * inputs["mu_par"]}


def symbolic():
    import sympy as sp
    mu_par = sp.Symbol("mu_par", positive=True)
    return {"mu_par_bar": mu_par / 4}
