"""Physical time corresponding to a count of fast live-imaging frames: N_frames * dt_frame."""


def run(inputs, settings):
    return {"t_dur": inputs["N_frames"] * inputs["dt_frame"]}


def symbolic():
    import sympy as sp

    N_frames, dt_frame = sp.symbols("N_frames dt_frame", positive=True)
    return {"t_dur": N_frames * dt_frame}
