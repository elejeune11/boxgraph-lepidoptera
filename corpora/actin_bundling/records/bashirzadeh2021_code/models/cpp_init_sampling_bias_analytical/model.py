"""Exact analytical moments of the initial-condition sampling distributions in:
crosslinker_sorting_in_GUVs/aster_initial_condition_grid_conda_freeze.py:
- Lines 47, 111-112: xkick, ykick ~ Uniform(-r_rand, +r_rand) -> std = r_rand / sqrt(3)
- Lines 156-158, 181-183: x ~ Uniform(-r_c, r_c), y|x ~ Uniform(-sqrt(r_c^2 - x^2), +sqrt(r_c^2 - x^2))
  Joint PDF: p(x, y) = 1 / (4 * r_c * sqrt(r_c^2 - x^2)) on x^2 + y^2 <= r_c^2.
  Exact moments:
    E[x^2] = r_c^2 / 3
    E[y^2] = E[(r_c^2 - x^2)/3] = (r_c^2 - r_c^2/3)/3 = 2*r_c^2 / 9
    E[x^2] / E[y^2] = 3/2 = 1.5 (vs 1.0 for any isotropic uniform disk or box)
    Density at (0, 0): N * p(0, 0) = (rho_xl * pi * r_c^2) / (4 * r_c^2) = (pi / 4) * rho_xl
    Fraction inside [-r_c/2, r_c/2]^2:
      (1 / (2*r_c)) * integral_{-r_c/2}^{r_c/2} (r_c / (2 * sqrt(r_c^2 - x^2))) dx
      = (1/2) * arcsin(1/2) = pi / 12.
"""

import math


def run(inputs, settings):
    r_rand = float(inputs["r_rand"])
    r_c = float(inputs["r_c"])
    rho_xl = float(inputs["rho_xl"])

    center_std_code = r_rand / math.sqrt(3.0)
    ex2_code = (r_c ** 2) / 3.0
    ey2_code = 2.0 * (r_c ** 2) / 9.0
    ex2_uniform_disk = (r_c ** 2) / 4.0
    anisotropy_code = ex2_code / ey2_code
    density_at_origin_code = (math.pi / 4.0) * rho_xl
    frac_in_rc_box_code = math.pi / 12.0

    return {
        "center_std_code": center_std_code,
        "ex2_code": ex2_code,
        "ey2_code": ey2_code,
        "ex2_uniform_disk": ex2_uniform_disk,
        "anisotropy_code": anisotropy_code,
        "density_at_origin_code": density_at_origin_code,
        "frac_in_rc_box_code": frac_in_rc_box_code,
    }


def symbolic():
    import sympy as sp

    r_rand, r_c, rho_xl = sp.symbols("r_rand r_c rho_xl", positive=True)
    return {
        "center_std_code": r_rand / sp.sqrt(3),
        "ex2_code": r_c ** 2 / 3,
        "ey2_code": sp.Rational(2, 9) * r_c ** 2,
        "ex2_uniform_disk": r_c ** 2 / 4,
        "anisotropy_code": sp.Rational(3, 2),
        "density_at_origin_code": (sp.pi / 4) * rho_xl,
        "frac_in_rc_box_code": sp.pi / 12,
    }
