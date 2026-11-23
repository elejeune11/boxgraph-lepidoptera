"""Motor ABD sliding velocity and filament traversal time from Eq. S4."""


def run(inputs, settings):
    v_walk = inputs["f_walk"] / inputs["zeta"]
    dalpha_dt = v_walk / inputs["l_f"]
    tau_walk = inputs["l_f"] / v_walk
    return {"v_walk": v_walk, "dalpha_dt": dalpha_dt, "tau_walk": tau_walk}


def symbolic():
    import sympy as sp

    f_walk, zeta, l_f = sp.symbols("f_walk zeta l_f", positive=True)
    v_walk = f_walk / zeta
    return {
        "v_walk": v_walk,
        "dalpha_dt": f_walk / (zeta * l_f),
        "tau_walk": l_f * zeta / f_walk,
    }
