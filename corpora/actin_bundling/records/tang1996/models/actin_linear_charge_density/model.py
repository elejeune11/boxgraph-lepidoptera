def run(inputs, settings):
    return {"sigma": inputs["q_net"] * inputs["rho_mon"]}


def symbolic():
    import sympy as sp

    q_net, rho_mon = sp.symbols("q_net rho_mon", positive=True)
    return {"sigma": q_net * rho_mon}
