import math


def run(inputs, settings):
    rho = float(inputs["rho"])
    tau = float(inputs["tau"])
    v = float(inputs["v"])
    T = float(inputs["T"])
    D_eff = tau * v**2 + T
    return {
        "M_cg": rho * D_eff,
        "mu_qs": math.log(rho) + 0.5 * math.log(D_eff),
    }


def symbolic():
    import sympy as sp

    rho, tau, v, T = sp.symbols("rho tau v T", positive=True)
    D_eff = tau * v**2 + T
    return {
        "M_cg": rho * D_eff,
        "mu_qs": sp.log(rho) + sp.Rational(1, 2) * sp.log(D_eff),
    }
