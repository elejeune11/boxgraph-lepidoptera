def run(inputs, settings):
    # Both W_bar and W_bar_c are negative attractive energies per length (W_bar < W_bar_c < 0 when bound).
    W_bar = inputs["W_bar"]
    W_bar_c = inputs["W_bar_c"]
    J = W_bar_c - W_bar if abs(W_bar) > abs(W_bar_c) else 0.0
    return {"J": max(0.0, J)}


def symbolic():
    import sympy as sp
    W_bar, W_bar_c = sp.symbols("W_bar W_bar_c", negative=True)
    return {"J": sp.Piecewise((W_bar_c - W_bar, -W_bar > -W_bar_c), (0, True))}
