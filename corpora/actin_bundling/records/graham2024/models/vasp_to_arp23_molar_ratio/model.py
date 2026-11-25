def run(inputs, settings):
    return {"ratio": inputs["C_vasp"] / inputs["C_arp"]}


def symbolic():
    import sympy as sp

    C_vasp, C_arp = sp.symbols("C_vasp C_arp", positive=True)
    return {"ratio": C_vasp / C_arp}
