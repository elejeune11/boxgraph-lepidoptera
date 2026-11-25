def run(inputs, settings):
    T_span = (inputs["N_frames"] - 1.0) * inputs["dt"]
    T_total = inputs["N_frames"] * inputs["dt"]
    return {"T_span": T_span, "T_total": T_total}


def symbolic():
    import sympy as sp

    N_frames, dt = sp.symbols("N_frames dt", positive=True)
    return {"T_span": (N_frames - 1) * dt, "T_total": N_frames * dt}
