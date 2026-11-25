import math

def run(inputs, settings):
    c0 = inputs["c_bar_0"]
    b = inputs["b_bar"]
    if c0 >= 0.0:
        S0 = 0.0
    else:
        S0 = math.sqrt(-c0 / b)
    return {"S_0": S0, "q_0_pos": 0.5 * S0, "q_0_neg": -0.5 * S0}
