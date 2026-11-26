def run(inputs, settings):
    # 1 revolution per minute = 1/60 revolution per second = 1/60 Hz
    return {"f_rot": inputs["omega_rpm"] / 60.0}


def symbolic():
    import sympy as sp

    omega_rpm = sp.symbols("omega_rpm", positive=True)
    return {"f_rot": omega_rpm / 60}
