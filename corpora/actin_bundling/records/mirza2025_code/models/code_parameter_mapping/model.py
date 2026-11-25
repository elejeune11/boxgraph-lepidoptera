"""Faithful Python port of the parameter-to-PDE mapping in mirza2025_code.

Ports:
- src/CortexNematic2D.cpp lines 225-258 (dparam packing: sus20/2, sus40/8, frank, etc.)
- src/AuxCortexNematic2D.cpp line 160 (q_int = -sqrt(-(4*sus20 + lambda_rot*heqb)/(32*sus40)))
- src/AuxCortexNematic2D.cpp line 420 (DqDq = 0.5*(|grad q1|^2 + |grad q2|^2) => L_bar = frank/2)
- src/AuxCortexNematic2D.cpp line 524 (-visc*(heqb - h) => k_bar_d = visc)
- src/AuxCortexNematic2D.cpp lines 527-528, 561-565 (3-parameter active stress => kappa)
- src/AuxCortexNematic2D.cpp lines 1070-1087, 1123-1140 (c_bar_0 and hidden log-barrier shift)
"""

import math


def run(inputs, settings):
    visc = float(inputs["visc"])
    cvisc = float(inputs["cvisc"])
    rvisc = float(inputs["rvisc"])
    frank = float(inputs["frank"])
    sus20_cfg = float(inputs["sus20"])
    sus40_cfg = float(inputs["sus40"])
    lambda_rot = float(inputs["lambda_rot"])
    lambda_iso = float(inputs["lambda_iso"])
    lambda_aniso = float(inputs["lambda_aniso"])
    crossLinker = float(inputs["crossLinker"])
    heqb = float(inputs["heqb"])
    k_lim = float(inputs["k_lim"])
    S2_max = float(inputs["S2_max"])
    r_out = float(inputs["r_out"])
    nx = float(inputs["nx"])

    # CortexNematic2D.cpp:226-227
    dparam_sus20 = sus20_cfg / 2.0
    dparam_sus40 = sus40_cfg / 8.0

    # Effective paper dimensionless groups from strong-form reduction of LS_CortexNematic_2D
    k_bar_d = visc
    L_bar_eff = 0.5 * frank
    lambda_o_bar = -lambda_rot
    c0_bar_bare = 4.0 * dparam_sus20 + lambda_rot * heqb  # = 2*sus20_cfg - lambda_o_bar*heqb
    c0_bar_with_barrier = c0_bar_bare + 4.0 * k_lim / S2_max

    kappa_eff = (lambda_aniso - crossLinker) / lambda_iso
    iso_S_excess = 0.5 * (lambda_aniso + crossLinker) / lambda_iso - 1.0

    # AuxCortexNematic2D.cpp:160: q_int = -sqrt(-(4*sus20 + lambda_rot*heqb)/(32*sus40))
    # (returns NaN in C++ when 4*sus20 + lambda_rot*heqb > 0, i.e., isotropic uniform state q_int = 0)
    radicand = -(4.0 * dparam_sus20 + lambda_rot * heqb) / (32.0 * dparam_sus40)
    q_int = -math.sqrt(radicand) if radicand > 0.0 else 0.0
    S_0 = 2.0 * abs(q_int)

    dx_elem = r_out / nx
    entropy_margin_2d = 2.0 * rvisc - cvisc * cvisc
    beta_bar_min_2d = -math.sqrt(2.0 * rvisc)

    return {
        "k_bar_d": k_bar_d,
        "L_bar_eff": L_bar_eff,
        "lambda_o_bar": lambda_o_bar,
        "c0_bar_bare": c0_bar_bare,
        "c0_bar_with_barrier": c0_bar_with_barrier,
        "kappa_eff": kappa_eff,
        "iso_S_excess": iso_S_excess,
        "q_int": q_int,
        "S_0": S_0,
        "dx_elem": dx_elem,
        "entropy_margin_2d": entropy_margin_2d,
        "beta_bar_min_2d": beta_bar_min_2d,
    }


def symbolic():
    import sympy as sp

    (
        visc,
        cvisc,
        rvisc,
        frank,
        sus20,
        sus40,
        lambda_rot,
        lambda_iso,
        lambda_aniso,
        crossLinker,
        heqb,
        k_lim,
        S2_max,
        r_out,
        nx,
    ) = sp.symbols(
        "visc cvisc rvisc frank sus20 sus40 lambda_rot lambda_iso lambda_aniso "
        "crossLinker heqb k_lim S2_max r_out nx",
        real=True,
    )
    c0_bare = 2 * sus20 + lambda_rot * heqb
    q_int_expr = sp.Piecewise((-sp.sqrt(-c0_bare / (4 * sus40)), c0_bare < 0), (0, True))
    return {
        "k_bar_d": visc,
        "L_bar_eff": frank / 2,
        "lambda_o_bar": -lambda_rot,
        "c0_bar_bare": c0_bare,
        "c0_bar_with_barrier": c0_bare + 4 * k_lim / S2_max,
        "kappa_eff": (lambda_aniso - crossLinker) / lambda_iso,
        "iso_S_excess": (lambda_aniso + crossLinker) / (2 * lambda_iso) - 1,
        "q_int": q_int_expr,
        "S_0": 2 * sp.Abs(q_int_expr),
        "dx_elem": r_out / nx,
        "entropy_margin_2d": 2 * rvisc - cvisc**2,
        "beta_bar_min_2d": -sp.sqrt(2 * rvisc),
    }
