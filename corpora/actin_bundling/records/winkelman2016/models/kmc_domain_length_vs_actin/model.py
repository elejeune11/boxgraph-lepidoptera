def run(inputs, settings):
    # Winkelman et al. (2016, Fig. S3B caption, p. 18 / SI-6) explicitly state that actin
    # filament growth is "not part of the model", so the cooperative kinetic lattice model
    # predicts a constant domain length independent of actin monomer concentration c_actin.
    _ = float(inputs["c_actin"])
    return {"l_dom": float(inputs["l_dom_const"])}


def symbolic():
    import sympy as sp

    c_actin, l_dom_const = sp.symbols("c_actin l_dom_const", positive=True)
    return {"l_dom": l_dom_const + 0 * c_actin}
