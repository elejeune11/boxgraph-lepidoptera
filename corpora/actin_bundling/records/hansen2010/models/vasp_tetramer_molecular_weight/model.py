def run(inputs, settings):
    return {"MW_tetramer": 4.0 * inputs["MW_monomer"]}


def symbolic():
    import sympy as sp

    MW_monomer = sp.symbols("MW_monomer", positive=True)
    return {"MW_tetramer": 4 * MW_monomer}
