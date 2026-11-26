import math


def _exact_moments(a, b, D_f, N_min):
    """Exact moments of P(N) propto exp(a*sqrt(N) - b*N) on [N_min, infinity).

    Substituting u = sqrt(N), dN = 2*u*du, s = sqrt(b), mu = a/(2*b), and
    t = s*(u - mu) reduces every moment to elementary + erfc integrals.
    """
    s = math.sqrt(b)
    inv_s = 1.0 / s
    mu = a / (2.0 * b)
    u0 = math.sqrt(max(0.0, N_min))
    t0 = s * (u0 - mu)
    e0 = math.exp(-t0 * t0)
    J0 = 0.5 * math.sqrt(math.pi) * math.erfc(t0)
    J1 = 0.5 * e0
    J2 = 0.5 * t0 * e0 + 0.5 * J0
    J3 = 0.5 * (t0 * t0 + 1.0) * e0

    I1 = (mu * J0 + inv_s * J1) * inv_s
    I2 = (mu**2 * J0 + 2.0 * mu * inv_s * J1 + inv_s**2 * J2) * inv_s
    I3 = (
        mu**3 * J0
        + 3.0 * mu**2 * inv_s * J1
        + 3.0 * mu * inv_s**2 * J2
        + inv_s**3 * J3
    ) * inv_s

    mean_N = I3 / I1
    mean_D = D_f * (I2 / I1)
    return mean_N, mean_D


def run(inputs, settings):
    c_on = inputs["two_kon_over_D"]
    c_off = inputs["koff_over_D"]
    D_f = float(inputs["D_f"])
    N_min = float(inputs["N_min"])

    xs = [float(x) for x in c_on["x"]]
    mean_Ns = []
    mean_Ds = []
    for a, b in zip(c_on["y"], c_off["y"]):
        mN, mD = _exact_moments(float(a), float(b), D_f, N_min)
        mean_Ns.append(mN)
        mean_Ds.append(mD)
    return {
        "N_b_mean": {"x": xs, "y": mean_Ns},
        "D_b_mean": {"x": xs, "y": mean_Ds},
    }
