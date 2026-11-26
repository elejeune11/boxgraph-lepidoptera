def run(inputs, settings):
    return {"N_walk": inputs["N_minifil"] * inputs["k_walk_cyl"] * inputs["t_min"]}

def symbolic():
    import sympy as sp
    N_minifil, k_walk_cyl, t_min = sp.symbols("N_minifil k_walk_cyl t_min", positive=True)
    return {"N_walk": N_minifil * k_walk_cyl * t_min}
