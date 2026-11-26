def run(inputs, settings):
    return {"alpha_prime": inputs["K"] / inputs["eps_b"]}


def symbolic():
    import sympy as sp

    K, eps_b = sp.symbols("K eps_b", positive=True)
    return {"alpha_prime": K / eps_b}
