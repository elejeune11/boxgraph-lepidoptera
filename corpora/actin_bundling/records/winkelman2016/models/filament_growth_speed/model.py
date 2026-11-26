def run(inputs, settings):
    return {"v": inputs["k_elong"] * inputs["d_sub"]}


def symbolic():
    import sympy as sp

    k_elong, d_sub = sp.symbols("k_elong d_sub", positive=True)
    return {"v": k_elong * d_sub}
