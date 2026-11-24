def run(inputs, settings):
    val = (
        inputs["kappa"]
        * inputs["Lambda"]
        * inputs["H0"]
        / inputs["ns"]
        * inputs["n"]
        * inputs["grad_H"]
    )
    return {"J_curv": val}


def symbolic():
    import sympy as sp

    kappa, Lambda, H0, ns, n, grad_H = sp.symbols(
        "kappa Lambda H0 ns n grad_H", real=True
    )
    return {"J_curv": (kappa * Lambda * H0 / ns) * n * grad_H}
