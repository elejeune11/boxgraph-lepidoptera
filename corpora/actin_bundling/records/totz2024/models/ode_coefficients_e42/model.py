def run(inputs, settings):
    nu = inputs["nu"]
    L = inputs["L"]
    h = inputs["h"]
    g = inputs["g"]
    c = inputs["c"]
    Y = inputs["Y"]
    p = inputs["p"]
    ratio = L / (2.0 * h)
    # Printed literally in SI Eq. E42 (evaluated in SI base units: m, Pa):
    # a = 12*(1 - nu^2)*((L/(2*h))^2 * g - c*Y*h^2)  [dimensionally inconsistent second term]
    a_ode = 12.0 * (1.0 - nu ** 2) * ((ratio ** 2) * g - c * Y * (h ** 2))
    b_ode = 12.0 * (1.0 - nu ** 2) * (p / Y) * (ratio ** 4)
    return {"a_ode": a_ode, "b_ode": b_ode}


def symbolic():
    import sympy as sp
    nu, L, h, g, c, Y, p = sp.symbols("nu L h g c Y p", positive=True)
    ratio = L / (2 * h)
    return {
        "a_ode": 12 * (1 - nu**2) * (ratio**2 * g - c * Y * h**2),
        "b_ode": 12 * (1 - nu**2) * (p / Y) * ratio**4,
    }
