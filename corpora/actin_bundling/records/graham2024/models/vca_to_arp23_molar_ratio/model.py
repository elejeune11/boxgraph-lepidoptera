def run(inputs, settings):
    return {"ratio": inputs["C_vca"] / inputs["C_arp"]}


def symbolic():
    import sympy as sp

    C_vca, C_arp = sp.symbols("C_vca C_arp", positive=True)
    return {"ratio": C_vca / C_arp}
