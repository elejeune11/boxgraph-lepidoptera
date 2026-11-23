import math


def run(inputs, settings):
    # Analytical integral of printed Phi(x, y) = [sinh((x-y-sign(x-y)l)/L) - sinh((x+y)/L)] / (2 sinh(l/L))
    # over y in [-l/2, l/2] for uniform density rho = rho_0:
    # int Phi dy = -L * sinh(x/L) / cosh(l/(2L))
    # Note: sigma_a [Pa*m] / eta [Pa*s*m] has unit [1/s]; times L [um] gives [um/s].
    L = inputs["L"]
    l = inputs["l"]
    x = inputs["x"]
    # Avoid overflow for large l/(2L): sinh(x/L)/cosh(l/(2L)) = (exp((x - l/2)/L) - exp((-x - l/2)/L)) / (1 + exp(-l/L))
    ratio = (math.exp((x - 0.5 * l) / L) - math.exp((-x - 0.5 * l) / L)) / (1.0 + math.exp(-l / L))
    v = -(inputs["sigma_a"] / inputs["eta"]) * L * ratio
    return {"v_unif": v}


def symbolic():
    import sympy as sp

    sigma_a, eta, L, l, x = sp.symbols("sigma_a eta L l x", positive=True)
    return {"v_unif": -(sigma_a / eta) * L * sp.sinh(x / L) / sp.cosh(l / (2 * L))}
