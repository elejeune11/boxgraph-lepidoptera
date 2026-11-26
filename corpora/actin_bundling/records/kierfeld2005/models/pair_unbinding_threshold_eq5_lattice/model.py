def run(inputs, settings):
    # For two identical filaments of persistence length L_p = 2*kappa/T, the reduced bending rigidity
    # is kappa_12^r = kappa / 2, so the reduced persistence length in Eq. (5) is L_p_rel = L_p / 2.
    L_p_rel = 0.5 * inputs["L_p"]
    l_a = inputs["l_a"]
    return {"W_bar_c": (1.0 / L_p_rel) * (L_p_rel / l_a) ** (2.0 / 3.0)}


def symbolic():
    import sympy as sp
    L_p, l_a = sp.symbols("L_p l_a", positive=True)
    L_p_rel = L_p / 2
    return {"W_bar_c": (1 / L_p_rel) * (L_p_rel / l_a) ** sp.Rational(2, 3)}
