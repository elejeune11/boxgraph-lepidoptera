def run(inputs, settings):
    # Verbatim printed equation on p. 2 col 2 has W_bar on both sides:
    # exp(-a_par * W_bar / T) = (1 - X_1) + X_1 * exp(-a_par * W_bar / T)
    # Rearranging gives (1 - X_1) * (exp(-a_par * W_bar / T) - 1) = 0,
    # whose unique real solution for X_1 < 1 is W_bar = 0 identically (independent of W).
    if inputs["X_1"] == 1.0:
        raise ValueError("Printed equation is indeterminate (0 = 0) at X_1 = 1")
    return {"W_bar": 0.0}


def symbolic():
    import sympy as sp
    return {"W_bar": sp.Integer(0)}
