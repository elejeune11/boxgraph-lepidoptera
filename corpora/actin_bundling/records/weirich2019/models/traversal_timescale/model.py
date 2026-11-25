def run(inputs, settings):
    # distance [um] * 1000 [nm/um] / v [nm/s] / 60 [s/min]
    return {"tau": (inputs["distance"] * 1000.0) / inputs["v"] / 60.0}


def symbolic():
    import sympy as sp
    distance, v = sp.symbols("distance v", positive=True)
    return {"tau": 1000 * distance / (60 * v)}
