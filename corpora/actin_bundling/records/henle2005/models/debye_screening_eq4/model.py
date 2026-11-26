import math

# 1 mol/L = N_A * 1e-24 nm^-3 (CODATA 2018 exact N_A = 6.02214076e23 mol^-1)
M_TO_NM3 = 0.602214076


def run(inputs, settings):
    l_B = inputs["l_B"]
    n_1 = inputs["n_1"] * M_TO_NM3
    n_q = inputs["n_q"] * M_TO_NM3
    q = inputs["q"]
    # Bulk electroneutrality: n_{-1} = n_1 + q * n_q => sum_s s^2 n_s = 2*n_1 + q*(q+1)*n_q
    ionic_sum = 2.0 * n_1 + q * (q + 1.0) * n_q
    return {"kappa": math.sqrt(4.0 * math.pi * l_B * ionic_sum)}


def symbolic():
    import sympy as sp

    l_B, n_1, n_q, q = sp.symbols("l_B n_1 n_q q", positive=True)
    ionic_sum = M_TO_NM3 * (2 * n_1 + q * (q + 1) * n_q)
    return {"kappa": sp.sqrt(4 * sp.pi * l_B * ionic_sum)}
