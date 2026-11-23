def run(inputs, settings):
    return {
        "R_0": (inputs["gamma"] ** 4)
        / (inputs["mu_perp_bar"] * (inputs["Sigma"] ** 3))
    }


def symbolic():
    import sympy as sp
    gamma, mu_perp_bar, Sigma = sp.symbols("gamma mu_perp_bar Sigma", positive=True)
    return {"R_0": gamma**4 / (mu_perp_bar * Sigma**3)}
