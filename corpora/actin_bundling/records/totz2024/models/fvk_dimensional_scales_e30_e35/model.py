def run(inputs, settings):
    L = inputs["L"]
    h = inputs["h"]
    Y = inputs["Y"]
    half_L = L / 2.0
    ratio = h / half_L
    return {
        "x_scale": half_L,
        "biharm_scale": 1.0 / (half_L ** 4),
        "g_scale": ratio ** 2,
        "w_scale": h,
        "chi_scale": Y * (h ** 2),
        "p_scale": (ratio ** 4) * Y,
    }


def symbolic():
    import sympy as sp
    L, h, Y = sp.symbols("L h Y", positive=True)
    half_L = L / 2
    return {
        "x_scale": half_L,
        "biharm_scale": 1 / half_L**4,
        "g_scale": (h / half_L)**2,
        "w_scale": h,
        "chi_scale": Y * h**2,
        "p_scale": (h / half_L)**4 * Y,
    }
