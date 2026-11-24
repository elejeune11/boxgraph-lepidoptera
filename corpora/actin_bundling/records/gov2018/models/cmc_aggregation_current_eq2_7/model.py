def run(inputs, settings):
    val = (
        2.0
        * inputs["J"]
        * inputs["Lambda"]
        / inputs["ns"]
        * inputs["n"]
        * (inputs["grad_n"] + inputs["grad3_n"] / inputs["ns"])
    )
    return {"J_agg": val}


def symbolic():
    import sympy as sp

    J, Lambda, ns, n, grad_n, grad3_n = sp.symbols(
        "J Lambda ns n grad_n grad3_n", real=True
    )
    return {"J_agg": (2 * J * Lambda / ns) * n * (grad_n + grad3_n / ns)}
