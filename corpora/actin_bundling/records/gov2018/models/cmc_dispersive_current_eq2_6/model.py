def run(inputs, settings):
    val = -(
        inputs["kappa"]
        * inputs["Lambda"]
        * (inputs["H0"] ** 2)
        / inputs["ns"]
        * inputs["n"]
        * inputs["grad_n"]
    )
    return {"J_disp": val}


def symbolic():
    import sympy as sp

    kappa, Lambda, H0, ns, n, grad_n = sp.symbols(
        "kappa Lambda H0 ns n grad_n", real=True
    )
    return {"J_disp": -(kappa * Lambda * H0 ** 2 / ns) * n * grad_n}
