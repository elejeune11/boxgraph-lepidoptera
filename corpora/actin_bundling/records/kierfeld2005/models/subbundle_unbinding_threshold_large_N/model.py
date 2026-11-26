def run(inputs, settings):
    return {
        "W_bar_c": inputs["W_bar_c_pair"]
        * (inputs["M"] ** (-1.0 / 3.0))
        / inputs["Delta_n_nn"]
    }


def symbolic():
    import sympy as sp
    W_bar_c_pair, M, Delta_n_nn = sp.symbols("W_bar_c_pair M Delta_n_nn", positive=True)
    return {"W_bar_c": W_bar_c_pair * M ** (-sp.Rational(1, 3)) / Delta_n_nn}
