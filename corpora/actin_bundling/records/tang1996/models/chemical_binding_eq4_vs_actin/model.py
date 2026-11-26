def run(inputs, settings):
    c_free = inputs["c_free"]
    v = inputs["V_P"] * inputs["q_net"] * inputs["c_actin"]
    return {"c_total": c_free + inputs["k"] * v * c_free}


def symbolic():
    import sympy as sp

    c_free, k, V_P, q_net, c_actin = sp.symbols("c_free k V_P q_net c_actin", positive=True)
    return {"c_total": c_free + k * V_P * q_net * c_actin * c_free}
