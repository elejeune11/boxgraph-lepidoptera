"""Initial membrane-bound VCA surface density assuming complete homogeneous adsorption inside a sphere."""

# 1 uM = 1e-6 mol/L = 1e-6 * N_A / 1e15 um^3 = 602.214076 molecules/um^3 (CODATA 2018 exact N_A)
NA_PER_UM3_PER_UM = 602.214076


def run(inputs, settings):
    c_vca = float(inputs["c_VCA"])  # uM
    r = float(inputs["R"])          # um
    # Sphere volume-to-surface-area ratio: ((4/3)*pi*R^3) / (4*pi*R^2) = R / 3
    v0 = (1.0 / 3.0) * c_vca * NA_PER_UM3_PER_UM * r
    return {"V0": v0}


def symbolic():
    import sympy as sp

    c_VCA, R = sp.symbols("c_VCA R", positive=True)
    return {"V0": sp.Rational(1, 3) * c_VCA * sp.Float("602.214076") * R}
