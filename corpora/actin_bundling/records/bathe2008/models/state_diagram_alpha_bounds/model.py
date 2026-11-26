"""Evaluate Eq. (1) coupling parameter alpha at [L_min, L_max] on the Fig. 5 box x-coordinates [N_min, N_max]."""


def run(inputs, settings):
    pref = inputs["k_x"] / (inputs["E_f_A_f"] * inputs["delta"])
    a_min = pref * (inputs["L_min"] ** 2)
    a_max = pref * (inputs["L_max"] ** 2)
    xs = list(inputs["box_ref"]["x"])
    return {"alpha_bounds": {"x": xs, "y": [a_min, a_max]}}
