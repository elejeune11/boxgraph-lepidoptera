def run(inputs, settings):
    r_grow_sites = inputs["mu_grow"] / inputs["l_b"]
    return {"rate_ratio": r_grow_sites / inputs["mu_on"]}


def symbolic():
    import sympy as sp

    mu_grow, l_b, mu_on = sp.symbols("mu_grow l_b mu_on", positive=True)
    return {"rate_ratio": mu_grow / (l_b * mu_on)}
