def run(inputs, settings):
    return {"frac_change": inputs["delta_q"] / inputs["q_net"]}


def symbolic():
    import sympy as sp

    delta_q, q_net = sp.symbols("delta_q q_net", positive=True)
    return {"frac_change": delta_q / q_net}
