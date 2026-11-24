"""Faithful Python port of CytoSim Hand kinetics and checks:
- src/sim/hand_prop.h:27 (POOL_UNATTACHED = 1)
- src/sim/hand_prop.cc:244-247, 296-301 (binding_prob, unbinding_rate_dt, unbinding_rate_log)
- src/sim/hand_prop.cc:314-333 (HandProp::checkStiffness for resting length len = 0)
- src/sim/hand_prop.cc:357-366 (HandProp::bindingSectionRate in 3D)
- src/sim/hand.h:286-302 (Hand::checkKramersDetachment)
- src/sim/couple.cc:175-178 (Couple::stepFF 1/2 coin flip)
"""
import math


def run(inputs, settings):
    k_bind = float(inputs["k_bind"])
    r_bind_um = float(inputs["binding_range"]) * 1e-3
    k_unbind0 = float(inputs["k_unbind0"])
    f_u = float(inputs["unbinding_force"])
    force = float(inputs["force"])
    dt = float(inputs["dt"])
    stiffness = float(inputs["stiffness"])
    kT = float(inputs["kT"])

    pool_unattached = 1.0  # src/sim/hand_prop.h:27
    br = k_bind * dt * pool_unattached
    binding_prob = -math.expm1(-br)  # src/sim/hand_prop.cc:247
    binding_prob_couple_first = 0.5 * binding_prob  # src/sim/couple.cc:175-178
    binding_section_rate_3d = math.pi * r_bind_um * r_bind_um * k_bind  # src/sim/hand_prop.cc:360

    unbinding_force_inv = 1.0 / f_u
    unbinding_rate_dt = k_unbind0 * dt
    unbinding_rate_log = math.log(unbinding_rate_dt)
    step_dec = math.exp(force * unbinding_force_inv + unbinding_rate_log)  # src/sim/hand.h:297-299
    k_unbind_loaded = step_dec / dt
    tau_cont = 1.0 / k_unbind_loaded
    # Exact discrete expectation E[ceil(Exp(1) / step_dec) * dt] = dt / (1 - exp(-step_dec))
    tau_disc = dt / (-math.expm1(-step_dec))

    # Port of HandProp::checkStiffness(stiff, len=0, ..., kT) (src/sim/hand_prop.cc:314-333)
    en_kT = (stiffness * r_bind_um * r_bind_um) / kT
    ap = math.exp(stiffness * r_bind_um * unbinding_force_inv)

    return {
        "binding_prob": binding_prob,
        "binding_prob_couple_first": binding_prob_couple_first,
        "binding_section_rate_3d": binding_section_rate_3d,
        "k_unbind_loaded": k_unbind_loaded,
        "tau_res_continuum": tau_cont,
        "tau_res_discrete": tau_disc,
        "check_stiffness_kT": en_kT,
        "check_stiffness_unbind_factor": ap,
    }


def symbolic():
    import sympy as sp

    k_bind, binding_range, k_unbind0, unbinding_force, force, dt, stiffness, kT = sp.symbols(
        "k_bind binding_range k_unbind0 unbinding_force force dt stiffness kT", positive=True
    )
    r_um = binding_range * sp.Rational(1, 1000)
    p_bind = 1 - sp.exp(-k_bind * dt)
    k_off = k_unbind0 * sp.exp(force / unbinding_force)
    return {
        "binding_prob": p_bind,
        "binding_prob_couple_first": sp.Rational(1, 2) * p_bind,
        "binding_section_rate_3d": sp.pi * r_um ** 2 * k_bind,
        "k_unbind_loaded": k_off,
        "tau_res_continuum": 1 / k_off,
        "tau_res_discrete": dt / (1 - sp.exp(-k_off * dt)),
        "check_stiffness_kT": stiffness * r_um ** 2 / kT,
        "check_stiffness_unbind_factor": sp.exp(stiffness * r_um / unbinding_force),
    }
