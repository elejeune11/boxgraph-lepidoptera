def run(inputs, settings):
    # Exact ratio from Eq. (5) using kappa_r_pair = kappa / 2 and kappa_r_M = M*(N-M)*kappa / N:
    # (kappa_r_pair / kappa_r_M)^(1/3) = (N / (2 * M * (N - M)))^(1/3).
    M = inputs["M"]
    N = inputs["N"]
    return {
        "W_bar_c": inputs["W_bar_c_pair"]
        * ((N / (2.0 * M * (N - M))) ** (1.0 / 3.0))
        / inputs["Delta_n_nn"]
    }


def symbolic():
    import sympy as sp
    W_bar_c_pair, M, N, Delta_n_nn = sp.symbols("W_bar_c_pair M N Delta_n_nn", positive=True)
    return {
        "W_bar_c": W_bar_c_pair * (N / (2 * M * (N - M))) ** sp.Rational(1, 3) / Delta_n_nn
    }
