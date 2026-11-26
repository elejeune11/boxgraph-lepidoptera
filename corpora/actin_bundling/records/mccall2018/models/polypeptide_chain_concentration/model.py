def run(inputs, settings):
    return {"c_chain": inputs["c_res"] / inputs["N_res"]}


def symbolic():
    import sympy as sp

    c_res, N_res = sp.symbols("c_res N_res", positive=True)
    return {"c_chain": c_res / N_res}
