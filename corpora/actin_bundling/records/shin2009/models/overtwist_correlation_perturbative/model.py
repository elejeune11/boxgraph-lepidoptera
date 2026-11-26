import cmath
import math

# Chirality-consistent right-handed mirror of the 6 printed 28/13 sites (m -> -m)
# so positive omega_0 = +12*pi/13 gives both V(l, l') > 0 and positive overtwist Im[ln g] > 0.
SITES_RH = ((0, 0), (4, -1), (9, 1), (14, 3), (19, -1), (24, 1))


def mean_overtwist_O2(beta_U: float, rho: float, xi_t: float, w0: float) -> float:
    z28 = cmath.exp(-28.0 / xi_t + 1j * 28.0 * w0)
    x = (beta_U * rho) ** 2
    ims = []
    for b in range(6):
        l0, m0 = SITES_RH[b]
        l1, m1 = SITES_RH[b + 1] if b < 5 else (28, 0)
        dl = l1 - l0

        # Exact periodic geometric sum S_j = sum_{j >= l1} exp(-(j-l1)/xi_t + i*(w0*(j-l1) - (phi_j - phi_1)))
        s_j_cell = 0j
        for k in range(6):
            ls, ms = SITES_RH[(b + 1 + k) % 6]
            d = (ls - l1) % 28
            dm = ms - m1
            s_j_cell += cmath.exp(-d / xi_t + 1j * (w0 * d - 2.0 * math.pi * dm / 6.0))
        S_j = s_j_cell / (1.0 - z28)

        # Exact periodic geometric sum S_k = sum_{k <= l0} exp(-(l0-k)/xi_t + i*(w0*(l0-k) - (phi_0 - phi_k)))
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
        ims.append(math.atan2(x * im_c, 1.0 + x * re_c))

    return sum(ims) / 6.0


def run(inputs, settings):
    beta_U = float(inputs["beta_U"])
    rho = float(inputs["rho"])
    xi_t = float(inputs["xi_t"])
    w0 = abs(float(inputs["omega_0"]))

    ov = mean_overtwist_O2(beta_U, rho, xi_t, w0)
    r1, r2 = rho * 0.99, rho * 1.01
    ov1 = mean_overtwist_O2(beta_U, r1, xi_t, w0)
    ov2 = mean_overtwist_O2(beta_U, r2, xi_t, w0)
    rho_exp = math.log(abs(ov2 / ov1)) / math.log(r2 / r1)
    return {"overtwist": ov, "rho_exponent": rho_exp}
