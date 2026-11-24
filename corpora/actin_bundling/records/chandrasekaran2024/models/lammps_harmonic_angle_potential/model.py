def run(inputs, settings):
    # LAMMPS angle_style harmonic convention: E = K * (theta - theta0)^2 with angles in radians
    return {"U_angle": inputs["K_angle"] * (inputs["theta"] - inputs["theta0"]) ** 2}


def symbolic():
    import sympy as sp

    K_angle, theta, theta0 = sp.symbols("K_angle theta theta0", positive=True)
    return {"U_angle": K_angle * (theta - theta0) ** 2}
