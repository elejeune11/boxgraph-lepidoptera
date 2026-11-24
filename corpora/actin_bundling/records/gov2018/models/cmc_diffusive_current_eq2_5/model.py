def run(inputs, settings):
    return {"J_diff": -inputs["D"] * inputs["grad_n"]}


def symbolic():
    import sympy as sp

    D, grad_n = sp.symbols("D grad_n", real=True)
    return {"J_diff": -D * grad_n}
