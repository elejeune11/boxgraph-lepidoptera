def run(inputs, settings):
    return {"d_linker": float(inputs["D"]) - float(inputs["d_actin"])}


def symbolic():
    import sympy as sp

    D, d_actin = sp.symbols("D d_actin", positive=True)
    return {"d_linker": D - d_actin}
