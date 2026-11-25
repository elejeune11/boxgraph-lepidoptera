def run(inputs, settings):
    return {"fold": inputs["P_vasp"] / inputs["P_arp"]}


def symbolic():
    import sympy as sp

    P_vasp, P_arp = sp.symbols("P_vasp P_arp", positive=True)
    return {"fold": P_vasp / P_arp}
