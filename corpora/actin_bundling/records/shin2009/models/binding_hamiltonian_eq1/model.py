def run(inputs, settings):
    eps_0 = float(inputs["eps_0"])
    k = float(inputs["k"])
    D = float(inputs["D"])
    a = float(inputs["a"])
    x = float(inputs["x"])
    # For monomer binding sites at R_i + a*S_i and R_j + a*S_j with R_j - R_i = D*D_hat_ij,
    # optimal alignment has S_i = +D_hat_ij and S_j = -D_hat_ij (x = 2), giving Delta_0 = D - 2*a.
    Delta_0 = D - 2.0 * a
    # Leading in-plane distance along lattice direction: Delta_ij = D - a*x + O(a^2/D)
    # => Delta_ij - Delta_0 = a*(2 - x), so (Delta_ij - Delta_0)^2 = a^2*(4 - 4*x + x^2)
    Delta_ij = D - a * x
    H_bond = -eps_0 + 0.5 * k * (Delta_ij - Delta_0) ** 2
    # Exact linear expansion coefficient of -a^2*x is 2*(D - Delta_0)/a = 4.0
    c_1_exact = 2.0 * (D - Delta_0) / a
    return {
        "Delta_0": Delta_0,
        "H_bond": H_bond,
        "c_1_exact": c_1_exact,
    }


def symbolic():
    import sympy as sp

    eps_0, k, D, a, x = sp.symbols("eps_0 k D a x", positive=True)
    Delta_0 = D - 2 * a
    Delta_ij = D - a * x
    return {
        "Delta_0": Delta_0,
        "H_bond": -eps_0 + sp.Rational(1, 2) * k * (Delta_ij - Delta_0) ** 2,
        "c_1_exact": sp.Integer(4),
    }
