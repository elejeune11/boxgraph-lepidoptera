import math


def run(inputs, settings):
    u = inputs["u"]
    Pe = inputs["Pe"]
    s0 = inputs["sigma0_hat"]
    # Numerically stable log(1 + exp(u))
    if u > 40.0:
        log1p_eu = u
    elif u < -40.0:
        log1p_eu = math.exp(u)
    else:
        log1p_eu = math.log1p(math.exp(u))
    V_dimless = -0.5 * u * u - Pe * s0 * u + Pe * (log1p_eu - math.log(2.0))
    # Numerically stable exp(u)/(1 + exp(u))^2 = 1 / (4 * cosh(u/2)^2)
    if abs(u) > 50.0:
        df_du = 0.0
    else:
        eu = math.exp(u)
        df_du = eu / ((1.0 + eu) ** 2)
    d2V = Pe * df_du - 1.0
    return {"V_dimless": V_dimless, "d2V_du2_dimless": d2V}


def symbolic():
    import sympy as sp

    u, sigma0_hat = sp.symbols("u sigma0_hat", real=True)
    Pe = sp.Symbol("Pe", positive=True)
    V_dimless = -u**2 / 2 - Pe * sigma0_hat * u + Pe * sp.log((1 + sp.exp(u)) / 2)
    d2V = Pe * sp.exp(u) / (1 + sp.exp(u)) ** 2 - 1
    return {"V_dimless": V_dimless, "d2V_du2_dimless": d2V}
