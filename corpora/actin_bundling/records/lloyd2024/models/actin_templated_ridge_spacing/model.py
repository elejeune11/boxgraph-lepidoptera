"""Adult cover scale ridge spacing from 1:1 F-actin bundle templating: s = W / N_actin (p. 3, SI p. 5)."""


def run(inputs, settings):
    return {"s": float(inputs["W"]) / float(inputs["N_actin"])}


def symbolic():
    import sympy as sp

    W, N_actin = sp.symbols("W N_actin", positive=True)
    return {"s": W / N_actin}
