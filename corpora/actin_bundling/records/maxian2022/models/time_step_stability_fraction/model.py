def run(inputs, settings):
    return {"frac_pct": 100.0 * inputs["dt"] / inputs["dt_max"]}


def symbolic():
    import sympy as sp

    dt, dt_max = sp.symbols("dt dt_max", positive=True)
    return {"frac_pct": 100 * dt / dt_max}
