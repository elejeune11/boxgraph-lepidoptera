def run(inputs, settings):
    total = inputs["dopc"] + inputs["dppc"] + inputs["chol"] + inputs["dotap"] + inputs["fitc"]
    return {"total_pct": total}


def symbolic():
    import sympy as sp

    dopc, dppc, chol, dotap, fitc = sp.symbols("dopc dppc chol dotap fitc", nonnegative=True)
    return {"total_pct": dopc + dppc + chol + dotap + fitc}
