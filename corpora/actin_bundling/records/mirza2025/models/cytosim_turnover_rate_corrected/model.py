def run(inputs, settings):
    return {"r_A": inputs["N_t"] / (inputs["N_s"] * inputs["dt"] * inputs["N_A"])}

def symbolic():
    import sympy as sp
    N_t, N_s, dt, N_A = sp.symbols("N_t N_s dt N_A", positive=True)
    return {"r_A": N_t / (N_s * dt * N_A)}
