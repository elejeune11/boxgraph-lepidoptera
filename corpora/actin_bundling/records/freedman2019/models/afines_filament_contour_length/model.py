def run(inputs, settings):
    return {"L": inputs["N_links"] * inputs["l_a"]}


def symbolic():
    import sympy as sp

    N_links, l_a = sp.symbols("N_links l_a", positive=True)
    return {"L": N_links * l_a}
