"""Map physical area bounds [A_min, A_max] (README.md line 5: [0.01, 0.1] um^2) to exact
interior integer pixel bounds [N_min, N_max] = [ceil(A_min / dx^2), floor(A_max / dx^2)]
corresponding to bwareafilt(BW_Image_1, [6 53]) on Nanocluster_motion.m line 40.
"""

import math


def run(inputs, settings):
    a_px = float(inputs["dx"]) ** 2
    n_min = math.ceil(float(inputs["A_min"]) / a_px)
    n_max = math.floor(float(inputs["A_max"]) / a_px)
    return {"N_min": float(n_min), "N_max": float(n_max)}
