"""Line-by-line Python port of AFINES KMC binding/unbinding and steric occlusion:
- src/motor.cpp:534-616 (motor::metropolis_prob)
- src/motor.cpp:621-639 (motor::alignment_penalty)
- src/motor.cpp:646-700 (motor::try_attach, motor::allowed_bind)
- src/motor.cpp:726-742 (motor::try_detach)
- src/filament.cpp:590-614 (filament::closest_attached_distance)
"""

import math


def _alignment_penalty(kalign, par_flag, theta):
    if kalign == 0.0:
        return 0.0
    c = max(-1.0, min(1.0, math.cos(theta)))
    if par_flag == 1:
        return kalign * (1.0 - c)
    if par_flag == -1:
        return kalign * (1.0 + c)
    return kalign * (1.0 - c) if c > 0.0 else kalign * (1.0 + c)


def run(inputs, settings):
    k_on = float(inputs["k_on"])
    k_off = float(inputs["k_off"])
    dt = float(inputs["dt"])
    temperature = float(inputs["temperature"])
    mk = float(inputs["mk"])
    mld = float(inputs["mld"])
    len_old = float(inputs["len_old"])
    len_new = float(inputs["len_new"])
    kalign = float(inputs["kalign"])
    par_flag = int(round(float(inputs["par_flag"])))
    theta = float(inputs["theta"])
    kb_motor = float(inputs["kb_motor"])
    same_filament = bool(int(round(float(inputs["same_filament"]))))
    same_spring = bool(int(round(float(inputs["same_spring"]))))
    occ = float(inputs["occ"])
    nearest_head_dist = float(inputs["nearest_head_dist"])

    # src/motor.cpp:537-539 stretching dE
    d_old = len_old - mld
    d_new = len_new - mld
    dE = 0.5 * mk * (d_new * d_new - d_old * d_old)

    # src/motor.cpp:580-608 alignment dE (assuming second head binding to doubly-bound state)
    dE += _alignment_penalty(kalign, par_flag, theta)

    # src/motor.cpp:615 Metropolis acceptance factor for attach (+dE) and reverse detach (-dE)
    p_metro_attach = 1.0 if dE <= 0.0 else math.exp(-dE / temperature)
    p_metro_detach = 1.0 if (-dE) <= 0.0 else math.exp(dE / temperature)

    # src/motor.cpp:696-700 motor::allowed_bind()
    # CRITICAL: when kb == 0 (a_m_bend = p_m_bend = 0 in sorting_kinetics.cfg),
    # both heads ARE allowed on the SAME filament as long as they are on different springs!
    if kb_motor > 0.0:
        allowed_topo = not same_filament
    else:
        allowed_topo = (not same_filament) or (not same_spring)

    # src/motor.cpp:684 steric occlusion gate
    allowed_occ = not (occ != 0.0 and nearest_head_dist < occ)

    p_attach_step = (k_on * dt) * p_metro_attach * (1.0 if (allowed_topo and allowed_occ) else 0.0)
    p_detach_step = (k_off * dt) * p_metro_detach
    K_effective = p_attach_step / p_detach_step if p_detach_step > 0.0 else 0.0

    return {
        "dE_attach": dE,
        "p_metro_attach": p_metro_attach,
        "p_metro_detach": p_metro_detach,
        "allowed_by_topology": 1.0 if allowed_topo else 0.0,
        "allowed_by_occ": 1.0 if allowed_occ else 0.0,
        "p_attach_step": p_attach_step,
        "p_detach_step": p_detach_step,
        "K_effective": K_effective,
    }
