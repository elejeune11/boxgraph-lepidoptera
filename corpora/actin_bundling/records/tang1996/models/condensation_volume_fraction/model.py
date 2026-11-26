def run(inputs, settings):
    return {"v": inputs["V_P"] * inputs["q_net"] * inputs["c_actin"]}


def symbolic():
    import sympy as sp

    V_P, q_net, c_actin = sp.symbols("V_P q_net c_actin", positive=True)
    return {"v": V_P * q_net * c_actin}
