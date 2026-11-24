def run(inputs, settings):
    return {"k_stretch": inputs["N_fil"] * inputs["k_stretch_fil"]}


def symbolic():
    import sympy as sp

    N_fil, k_stretch_fil = sp.symbols("N_fil k_stretch_fil", positive=True)
    return {"k_stretch": N_fil * k_stretch_fil}
