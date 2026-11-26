_R_KCAL_MOL_K = 8.31446261815324 / 4184.0


def run(inputs, settings):
    return {"E_molar": float(inputs["eps_perp"] * _R_KCAL_MOL_K * inputs["T"])}


def symbolic():
    import sympy as sp

    eps_perp, T = sp.symbols("eps_perp T", real=True)
    return {"E_molar": eps_perp * _R_KCAL_MOL_K * T}
