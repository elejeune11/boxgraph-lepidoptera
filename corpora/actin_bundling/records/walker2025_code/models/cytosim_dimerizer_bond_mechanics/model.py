"""Port of custom CytoSim Dimerizer double-link mechanics and discrete binding/unbinding."""

import math


def run(inputs, settings):
    k_single = float(inputs["k_single"])
    dx = float(inputs["dx"])
    k_bind = float(inputs["k_bind"])
    k_unbind_0 = float(inputs["k_unbind_0"])
    F_unbind = float(inputs["F_unbind"])
    dt = float(inputs["dt"])

    # Both attached dimerizer Wrists execute anchor.addLink(..., prop->stiffness) in simul_solve.cc:239
    k_eff = 2.0 * k_single
    F_mech = k_eff * dx
    # Wrist::force() computes prop->stiffness * d = k_single * dx (half of F_mech)
    F_bell = k_single * dx

    k_unbind_code = k_unbind_0 * math.exp(F_bell / F_unbind)
    k_unbind_at_Fmech = k_unbind_0 * math.exp(F_mech / F_unbind)

    p_bind_single = -math.expm1(-k_bind * dt)
    p_bind_pair = 1.0 - (1.0 - p_bind_single) ** 2

    return {
        "k_eff": k_eff,
        "F_mech": F_mech,
        "F_bell": F_bell,
        "k_unbind_code": k_unbind_code,
        "k_unbind_at_Fmech": k_unbind_at_Fmech,
        "p_bind_single": p_bind_single,
        "p_bind_pair": p_bind_pair,
    }


def symbolic():
    import sympy as sp

    k_single, dx, k_bind, k_unbind_0, F_unbind, dt = sp.symbols(
        "k_single dx k_bind k_unbind_0 F_unbind dt", positive=True
    )
    k_eff = 2 * k_single
    F_mech = k_eff * dx
    F_bell = k_single * dx
    k_unbind_code = k_unbind_0 * sp.exp(F_bell / F_unbind)
    k_unbind_at_Fmech = k_unbind_0 * sp.exp(F_mech / F_unbind)
    p_bind_single = 1 - sp.exp(-k_bind * dt)
    p_bind_pair = 1 - sp.exp(-2 * k_bind * dt)
    return {
        "k_eff": k_eff,
        "F_mech": F_mech,
        "F_bell": F_bell,
        "k_unbind_code": k_unbind_code,
        "k_unbind_at_Fmech": k_unbind_at_Fmech,
        "p_bind_single": p_bind_single,
        "p_bind_pair": p_bind_pair,
    }
