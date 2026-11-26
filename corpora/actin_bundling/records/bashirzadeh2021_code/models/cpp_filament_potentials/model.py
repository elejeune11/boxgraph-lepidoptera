"""Line-by-line Python port of AFINES filament stretching and bending potentials:
- src/spring.cpp:54-58, 185-188 (spring::update_force, spring::get_stretching_energy)
- prog/network.cpp:449 (link_bending_stiffness = polymer_bending_modulus / link_length)
- src/filament.cpp:320-345 (filament::update_bending)
- src/potentials.cpp:91-117 (bend_harmonic) with include/globals.h:61 (maxSmallAngle = 0.001)
"""

import math


def run(inputs, settings):
    k_s = float(inputs["k_s"])
    l_0 = float(inputs["l_0"])
    l_cur = float(inputs["l_cur"])
    kappa_b = float(inputs["kappa_b"])
    theta = float(inputs["theta"])
    max_small_angle = float(inputs["max_small_angle"])

    # src/spring.cpp:56, 187
    F_stretch = k_s * (l_cur - l_0)
    U_stretch = 0.5 * k_s * (l_cur - l_0) ** 2

    # prog/network.cpp:449 & src/potentials.cpp:91-117
    k_b = kappa_b / l_0
    c = math.cos(theta)
    if c > 1.0:
        c = 1.0
    if c < -1.0:
        c = -1.0
    s_raw = math.sqrt(max(0.0, 1.0 - c * c))
    s = max(s_raw, max_small_angle)
    ang = math.acos(c)
    U_bend = 0.5 * k_b * ang * ang
    # Restoring torque magnitude |delr2 x (a * tmp.force2)| = (k_b * ang / s) * s_raw
    tau_bend = (k_b * ang / s) * s_raw

    return {
        "U_stretch": U_stretch,
        "F_stretch": F_stretch,
        "k_b": k_b,
        "U_bend": U_bend,
        "tau_bend": tau_bend,
    }
