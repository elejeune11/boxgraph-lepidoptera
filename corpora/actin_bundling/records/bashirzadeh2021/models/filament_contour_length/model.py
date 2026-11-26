def run(inputs, settings):
    return {"L_f": (inputs["N_beads"] - 1.0) * inputs["l_link"]}


def symbolic():
    import sympy as sp

    N_beads, l_link = sp.symbols("N_beads l_link", positive=True)
    return {"L_f": (N_beads - 1) * l_link}
