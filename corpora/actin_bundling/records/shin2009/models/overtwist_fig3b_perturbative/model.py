import cmath
import math

SITES_RH = ((0, 0), (4, -1), (9, 1), (14, 3), (19, -1), (24, 1))


def mean_overtwist_O2_xsq(x_sq: float, xi_t: float, w0: float) -> float:
    z28 = cmath.exp(-28.0 / xi_t + 1j * 28.0 * w0)
    ims = []
    for b in range(6):
        l0, m0 = SITES_RH[b]
        l1, m1 = SITES_RH[b + 1] if b < 5 else (28, 0)
        dl = l1 - l0

        s_j_cell = 0j
        for k in range(6):
            ls, ms = SITES_RH[(b + 1 + k) % 6]
            d = (ls - l1) % 28
            dm = ms - m1
            s_j_cell += cmath.exp(-d / xi_t + 1j * (w0 * d - 2.0 * math.pi * dm / 6.0))
        S_j = s_j_cell / (1.0 - z28)

        s_k_cell = 0j
        for k in range(6):
            ls, ms = SITES_RH[(b - k) % 6]
            d = (l0 - ls) % 28
            dm = m0 - ms
            s_k_cell += cmath.exp(-d / xi_t + 1j * (w0 * d - 2.0 * math.pi * dm / 6.0))
        S_k = s_k_cell / (1.0 - z28)

        S = S_j * S_k * cmath.exp(-dl / xi_t + 1j * (w0 * dl - 2.0 * math.pi * (m1 - m0) / 6.0))
        re_c = 0.5 * (math.cosh(2.0 * dl / xi_t) - 1.0) * S.real
        im_c = -0.5 * math.sinh(2.0 * dl / xi_t) * S.imag
        ims.append(math.atan2(x_sq * im_c, 1.0 + x_sq * re_c))
    return sum(ims) / 6.0


def run(inputs, settings):
    x_sq = float(inputs["x_sq"])
    xi_t = float(inputs["xi_t"])
    w0 = abs(float(inputs["omega_0"]))
    return {"overtwist": mean_overtwist_O2_xsq(x_sq, xi_t, w0)}
