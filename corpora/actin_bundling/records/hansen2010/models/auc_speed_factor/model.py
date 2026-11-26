def run(inputs, settings):
    ratio = inputs["omega"] / inputs["omega_ref"]
    return {"speed_factor": ratio * ratio}


def symbolic():
    import sympy as sp

    omega, omega_ref = sp.symbols("omega omega_ref", positive=True)
    return {"speed_factor": (omega / omega_ref) ** 2}
