def run(inputs, settings):
    C_tot = float(inputs["C_tot"])
    r = float(inputs["r_Lpd_VASP"])
    C_VASP = C_tot / (1.0 + r)
    C_Lpd = C_tot * r / (1.0 + r)
    return {"C_VASP": C_VASP, "C_Lpd": C_Lpd}


def symbolic():
    import sympy as sp

    C_tot, r_Lpd_VASP = sp.symbols("C_tot r_Lpd_VASP", positive=True)
    return {
        "C_VASP": C_tot / (1 + r_Lpd_VASP),
        "C_Lpd": C_tot * r_Lpd_VASP / (1 + r_Lpd_VASP),
    }
