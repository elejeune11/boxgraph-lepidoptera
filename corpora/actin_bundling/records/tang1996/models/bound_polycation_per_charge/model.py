def run(inputs, settings):
    return {"theta_18": inputs["slope"] / inputs["q_net"]}


def symbolic():
    import sympy as sp

    slope, q_net = sp.symbols("slope q_net", positive=True)
    return {"theta_18": slope / q_net}
