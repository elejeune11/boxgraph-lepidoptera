def run(inputs, settings):
    return {"r_bind": inputs["r_Arp"] + inputs["d_actin"]}


def symbolic():
    import sympy as sp

    r_Arp, d_actin = sp.symbols("r_Arp d_actin", positive=True)
    return {"r_bind": r_Arp + d_actin}
