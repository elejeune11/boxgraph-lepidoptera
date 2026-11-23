def run(inputs, settings):
    # l_p in um -> multiply by 1000 to convert to nm
    return {
        "Omega": (inputs["gamma"] / (inputs["kBT"] * 1000.0 * inputs["l_p"]))
        ** (1.0 / 3.0)
    }


def symbolic():
    import sympy as sp
    gamma, kBT, l_p = sp.symbols("gamma kBT l_p", positive=True)
    return {"Omega": (gamma / (kBT * 1000 * l_p)) ** sp.Rational(1, 3)}
