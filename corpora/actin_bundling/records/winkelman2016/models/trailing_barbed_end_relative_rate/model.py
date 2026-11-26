def run(inputs, settings):
    v_lead = inputs["v_leading"]
    v_trail = inputs["v_trailing"]
    return {
        "f_rel": 100.0 * v_trail / v_lead,
        "f_red": 100.0 * (v_lead - v_trail) / v_lead,
    }


def symbolic():
    import sympy as sp

    v_leading, v_trailing = sp.symbols("v_leading v_trailing", positive=True)
    return {
        "f_rel": 100 * v_trailing / v_leading,
        "f_red": 100 * (v_leading - v_trailing) / v_leading,
    }
