def run(inputs, settings):
    rho = inputs["rho"]
    alpha = inputs["alpha"]
    return {
        "Delta": alpha**2 - rho**2,
        "Sigma": (rho + alpha) ** 2,
    }


def symbolic():
    import sympy as sp

    rho, alpha = sp.symbols("rho alpha", real=True)
    return {
        "Delta": alpha**2 - rho**2,
        "Sigma": (rho + alpha) ** 2,
    }
