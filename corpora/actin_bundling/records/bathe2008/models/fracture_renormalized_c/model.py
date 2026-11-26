"""Effective Eq. (2) mode factor c_eff propto (phi / phi*)^2 in the high-fracture-density regime (Fig. S2)."""


def run(inputs, settings):
    return {"c_eff": inputs["phi_ratio_sq"]}


def symbolic():
    import sympy as sp

    phi_ratio_sq = sp.symbols("phi_ratio_sq", positive=True)
    return {"c_eff": phi_ratio_sq}
