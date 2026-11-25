def run(inputs, settings):
    return {"N_dimer_max": float(inputs["N_monomer"]) / 2.0}


def symbolic():
    import sympy as sp

    N_monomer = sp.Symbol("N_monomer", positive=True)
    return {"N_dimer_max": N_monomer / 2}
