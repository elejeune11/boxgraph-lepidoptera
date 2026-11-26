def run(inputs, settings):
    return {"r_xl": inputs["C_xl"] / inputs["C_actin"]}


def symbolic():
    import sympy as sp

    C_xl, C_actin = sp.symbols("C_xl C_actin", positive=True)
    return {"r_xl": C_xl / C_actin}
