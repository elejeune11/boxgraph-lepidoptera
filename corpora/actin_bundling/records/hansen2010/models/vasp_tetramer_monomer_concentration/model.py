def run(inputs, settings):
    return {"C_tetramer": inputs["C_monomer"] / 4.0}


def symbolic():
    import sympy as sp

    C_monomer = sp.symbols("C_monomer", positive=True)
    return {"C_tetramer": C_monomer / 4}
