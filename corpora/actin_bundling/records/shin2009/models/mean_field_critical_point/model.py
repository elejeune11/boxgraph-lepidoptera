import math
from scipy.optimize import brentq


def dlnz_drho(rho: float, u: float) -> float:
    return 1.0 / rho + 1.0 / (1.0 - rho) - u


def d2lnz_drho2(rho: float) -> float:
    return -1.0 / (rho * rho) + 1.0 / ((1.0 - rho) ** 2)


def run(inputs, settings):
    tol = float(settings.get("tol", 1.0e-12))
    dr = float(settings.get("delta_rho", 1.0e-3))

    # Critical point: d(ln z)/drho = 0 and d^2(ln z)/drho^2 = 0
    rho_c = float(brentq(d2lnz_drho2, 0.05, 0.95, xtol=tol, rtol=tol))
    u_c = float(brentq(lambda u: dlnz_drho(rho_c, u), 1.0, 10.0, xtol=tol, rtol=tol))
    lnz_c = math.log(rho_c / (1.0 - rho_c)) - u_c * rho_c
    z_c = math.exp(lnz_c)

    # Evaluate susceptibility chi_rho = 1 / (z * d(ln z)/drho) at rho_c + dr and rho_c + 2*dr
    r1 = rho_c + dr
    r2 = rho_c + 2.0 * dr
    dlnz1 = math.log(r1 / (1.0 - r1)) - u_c * r1 - lnz_c
    dlnz2 = math.log(r2 / (1.0 - r2)) - u_c * r2 - lnz_c
    dz1 = abs(z_c * math.expm1(dlnz1))
    dz2 = abs(z_c * math.expm1(dlnz2))
    chi1 = 1.0 / (math.exp(lnz_c + dlnz1) * dlnz_drho(r1, u_c))
    chi2 = 1.0 / (math.exp(lnz_c + dlnz2) * dlnz_drho(r2, u_c))
    chi_exponent = math.log(chi2 / chi1) / math.log(dz2 / dz1)

    return {
        "u_c": u_c,
        "rho_c": rho_c,
        "z_c": z_c,
        "chi_exponent": chi_exponent,
    }
