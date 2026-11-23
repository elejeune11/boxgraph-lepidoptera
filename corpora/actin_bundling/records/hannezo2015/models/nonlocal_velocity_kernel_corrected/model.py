def run(inputs, settings):
    # With +sinh((x+y)/L), int_{-l/2}^{l/2} Phi_corrected(x, y) dy = 0 identically for all x, L, l
    return {"v_unif": 0.0}


def symbolic():
    import sympy as sp

    return {"v_unif": sp.Integer(0)}
