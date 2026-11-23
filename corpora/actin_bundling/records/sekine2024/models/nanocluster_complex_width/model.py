"""Widths of lifeact::GFP-labeled actin filament and alpha-actinin crosslinked dimer (Supp. Fig. 6a)."""


def run(inputs, settings):
    d_actin = inputs["d_actin"]
    d_gfp = inputs["d_GFP"]
    L_actinin = inputs["L_actinin"]
    return {
        "w_single": d_actin + 2.0 * d_gfp,
        "w_dimer": 2.0 * d_actin + 2.0 * d_gfp + L_actinin,
    }


def symbolic():
    import sympy as sp

    d_actin, d_GFP, L_actinin = sp.symbols("d_actin d_GFP L_actinin", positive=True)
    return {
        "w_single": d_actin + 2 * d_GFP,
        "w_dimer": 2 * d_actin + 2 * d_GFP + L_actinin,
    }
