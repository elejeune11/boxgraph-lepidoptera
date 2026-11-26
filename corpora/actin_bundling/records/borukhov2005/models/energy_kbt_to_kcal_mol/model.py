# Molar gas constant R = N_A * k_B in kcal / (mol * K) (thermochemical calorie = 4.184 J):
# R = 8.31446261815324 / 4184.0 = 1.9872042586408316e-3 kcal / (mol * K)
_R_KCAL_MOL_K = 8.31446261815324 / 4184.0


def run(inputs, settings):
    return {"E_molar": float(inputs["eps"] * _R_KCAL_MOL_K * inputs["T"])}


def symbolic():
    import sympy as sp

    eps, T = sp.symbols("eps T", real=True)
    return {"E_molar": eps * _R_KCAL_MOL_K * T}
