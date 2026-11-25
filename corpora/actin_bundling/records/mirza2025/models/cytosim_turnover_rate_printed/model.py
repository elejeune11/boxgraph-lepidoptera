def run(inputs, settings):
    return {"r_A_per_step": inputs["N_t"] / (inputs["N_s"] * inputs["N_A"])}

def symbolic():
    import sympy as sp
    N_t, N_s, N_A = sp.symbols("N_t N_s N_A", positive=True)
    return {"r_A_per_step": N_t / (N_s * N_A)}
