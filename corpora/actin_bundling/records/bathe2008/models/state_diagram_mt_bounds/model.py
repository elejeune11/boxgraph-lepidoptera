"""Evaluate alpha = G_shear * L^2 / (E_f * A_f) at [L_min, L_max] on Fig. 5g box x-coordinates."""


def run(inputs, settings):
    pref = inputs["G_shear"] / inputs["E_f_A_f"]
    a_min = pref * (inputs["L_min"] ** 2)
    a_max = pref * (inputs["L_max"] ** 2)
    xs = list(inputs["box_ref"]["x"])
    return {"alpha_bounds": {"x": xs, "y": [a_min, a_max]}}
