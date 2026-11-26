def run(inputs, settings):
    return {"a_mu": inputs["gamma_mu"] * inputs["c_mu"]}


def symbolic():
    import sympy as sp
    gamma_mu, c_mu = sp.symbols("gamma_mu c_mu", positive=True)
    return {"a_mu": gamma_mu * c_mu}
