def run(inputs, settings):
    return {"s": inputs["k"] * inputs["V_P"] * inputs["q_net"] * inputs["c_free"]}


def symbolic():
    import sympy as sp

    k, V_P, q_net, c_free = sp.symbols("k V_P q_net c_free", positive=True)
    return {"s": k * V_P * q_net * c_free}
