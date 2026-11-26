def run(inputs, settings):
    return {
        "v": inputs["v_0"]
        + inputs["delta_v_max"] * inputs["C_vasp"] / (inputs["K_d"] + inputs["C_vasp"])
    }


def symbolic():
    import sympy as sp

    v_0, delta_v_max, C_vasp, K_d = sp.symbols("v_0 delta_v_max C_vasp K_d", positive=True)
    return {"v": v_0 + delta_v_max * C_vasp / (K_d + C_vasp)}
