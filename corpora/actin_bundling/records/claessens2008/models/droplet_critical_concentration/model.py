N_A = 6.02214076e23  # 1/mol (exact SI)


def run(inputs, settings):
    # Inverse of n_f = (1/6) * D_d^2 * (c_a / M_actin) * N_A * d * 1e-18
    return {
        "c_a": 6.0 * inputs["n_f"] * inputs["M_actin"] / (N_A * 1e-18 * inputs["d"] * (inputs["D_d"] ** 2))
    }


def symbolic():
    import sympy as sp

    n_f, D_d, M_actin, d = sp.symbols("n_f D_d M_actin d", positive=True)
    coeff = sp.Float(N_A * 1e-18)
    return {"c_a": 6 * n_f * M_actin / (coeff * d * D_d**2)}
