"""Fiber coupling parameter alpha = (k_x / delta) * L^2 / (E_f * A_f)."""


def run(inputs, settings):
    # Pa * m^2 / N = (N/m^2) * m^2 / N = dimensionless
    return {"alpha": inputs["G_shear"] * (inputs["L"] ** 2) / inputs["E_f_A_f"]}


def symbolic():
    import sympy as sp

    G_shear, E_f_A_f, L = sp.symbols("G_shear E_f_A_f L", positive=True)
    return {"alpha": G_shear * L**2 / E_f_A_f}
