def run(inputs, settings):
    t = float(inputs["t"])
    chi_start = float(inputs["chi_start"])
    chi_end = float(inputs["chi_end"])
    t_ramp = float(inputs["t_ramp"])
    if t <= 0.0:
        chi = chi_start
    elif t >= t_ramp:
        chi = chi_end
    else:
        chi = chi_start + (chi_end - chi_start) * (t / t_ramp)
    return {"chi": chi}
