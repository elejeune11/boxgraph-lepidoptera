import math
from scipy.optimize import brentq


def run(inputs, settings):
    f_0 = inputs["f_0"]
    lam3 = inputs["lambda_3perp"]
    Sigma = inputs["Sigma"]
    rtol = float(settings.get("rtol", 1e-8))

    sigma_c = f_0 / ((2.0 ** (4.0 / 3.0)) * math.pi * (lam3 ** (1.0 / 3.0)))
    if Sigma > sigma_c * (1.0 + 1e-12):
        return {"R_0": None}
    if abs(Sigma - sigma_c) <= 1e-12 * sigma_c:
        return {"R_0": lam3}

    # Stationarity of g(u) = -u^(4/3)/(1+u^2)^(1/3) + s*u where u = R/lambda_3perp, s = 2*pi*Sigma*lambda_3perp^(1/3)/f_0
    # g'(u) = 0 <=> u^(1/3)*(4 + 2*u^2) / (3*(1+u^2)^(4/3)) = s, strictly decreasing on u in [1, oo)
    s = 2.0 * math.pi * Sigma * (lam3 ** (1.0 / 3.0)) / f_0

    def resid(u):
        return (u ** (1.0 / 3.0)) * (4.0 + 2.0 * u**2) / (3.0 * ((1.0 + u**2) ** (4.0 / 3.0))) - s

    u_lo = 1.0
    u_hi = max(2.0, (2.0 / (3.0 * s)) ** 3 * 2.0)
    while resid(u_hi) > 0.0:
        u_hi *= 2.0
    u_star = brentq(resid, u_lo, u_hi, rtol=rtol, xtol=1e-300)
    return {"R_0": u_star * lam3}
