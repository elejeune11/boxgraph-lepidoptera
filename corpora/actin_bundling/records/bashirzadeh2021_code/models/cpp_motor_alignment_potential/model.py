"""Line-by-line Python port of:
- src/potentials.cpp:38-89 (bend_angle, angle)
- src/motor.cpp:392-435 (motor::update_alignment)
- src/motor.cpp:621-639 (motor::alignment_penalty)
- src/motor.cpp:580-615 (alignment contribution to motor::metropolis_prob)
"""

import math


def run(inputs, settings):
    kalign = float(inputs["kalign"])
    par_flag = int(round(float(inputs["par_flag"])))
    theta = float(inputs["theta"])
    temperature = float(inputs["temperature"])

    # Unit vectors delr0 = (1, 0), delr1 = (cos(theta), sin(theta))
    dx0, dy0 = 1.0, 0.0
    dx1, dy1 = math.cos(theta), math.sin(theta)

    # bend_angle(delr0, delr1) from src/potentials.cpp:38-69
    rsq0 = dx0 * dx0 + dy0 * dy0
    r0 = math.sqrt(rsq0)
    rsq1 = dx1 * dx1 + dy1 * dy1
    r1 = math.sqrt(rsq1)
    c = (dx0 * dx1 + dy0 * dy1) / (r0 * r1)
    if c > 1.0:
        c = 1.0
    if c < -1.0:
        c = -1.0

    a11 = c / rsq0
    a12 = -1.0 / (r0 * r1)
    a22 = c / rsq1
    f1x = a11 * dx0 + a12 * dx1
    f1y = a11 * dy0 + a12 * dy1
    f2x = a22 * dx1 + a12 * dx0
    f2y = a22 * dy1 + a12 * dy0

    # motor::update_alignment() / motor::alignment_penalty() (src/motor.cpp:407-425, 628-638)
    if kalign == 0.0:
        a = 0.0
        align_eng = 0.0
    elif par_flag == 1:
        a = -kalign
        align_eng = kalign * (1.0 - c)
    elif par_flag == -1:
        a = kalign
        align_eng = kalign * (1.0 + c)
    else:
        if c > 0.0:
            a = -kalign
            align_eng = kalign * (1.0 - c)
        else:
            a = kalign
            align_eng = kalign * (1.0 + c)

    # Restoring torque magnitude |delr1 x (a * force2)| = |a| * |sin(theta)|
    fx_on_1 = a * f2x
    fy_on_1 = a * f2y
    tau_align = abs(dx1 * fy_on_1 - dy1 * fx_on_1)

    p_metro = 1.0 if align_eng <= 0.0 else math.exp(-align_eng / temperature)

    return {
        "U_align": align_eng,
        "tau_align": tau_align,
        "p_metro_align": p_metro,
    }
