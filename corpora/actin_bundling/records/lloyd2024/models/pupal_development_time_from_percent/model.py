"""Elapsed pupal development time t = (p_dev / 100) * T_pupal (p. 12)."""


def run(inputs, settings):
    # Port p_dev is in percent (e.g. 50 -> 50.0); T_pupal is in day.
    t_days = (float(inputs["p_dev"]) / 100.0) * float(inputs["T_pupal"])
    return {"t_days": t_days, "t_hours": t_days * 24.0}


def symbolic():
    import sympy as sp

    p_dev, T_pupal = sp.symbols("p_dev T_pupal", positive=True)
    t_days = (p_dev / 100) * T_pupal
    return {"t_days": t_days, "t_hours": 24 * t_days}
