def run(inputs, settings):
    return {"alpha": inputs["k_prime"] / inputs["eps_b"]}


def symbolic():
    import sympy as sp

    k_prime, eps_b = sp.symbols("k_prime eps_b", positive=True)
    return {"alpha": k_prime / eps_b}
