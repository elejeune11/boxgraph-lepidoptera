def run(inputs, settings):
    return {"c_res": inputs["N_res"] * inputs["c_chain"]}


def symbolic():
    import sympy as sp

    c_chain, N_res = sp.symbols("c_chain N_res", positive=True)
    return {"c_res": N_res * c_chain}
