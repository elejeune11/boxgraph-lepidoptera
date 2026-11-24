import math


def run(inputs, settings):
    rho = inputs["rho"]
    alpha = inputs["alpha"]
    Sigma = inputs["Sigma"]
    if Sigma < 0.0:
        return {"Delta": None}
    s = math.sqrt(Sigma)
    # Eq. (23): Delta(Sigma) = -2*rho*sqrt(Sigma) + Sigma for alpha > -rho,
    #           Delta(Sigma) = +2*rho*sqrt(Sigma) + Sigma for alpha < -rho.
    sign = -1.0 if alpha >= -rho else 1.0
    return {"Delta": sign * 2.0 * rho * s + Sigma}


def symbolic():
    import sympy as sp

    rho, alpha, Sigma = sp.symbols("rho alpha Sigma", real=True)
    sign = sp.Piecewise((-1, alpha >= -rho), (1, True))
    return {"Delta": sign * 2 * rho * sp.sqrt(Sigma) + Sigma}
