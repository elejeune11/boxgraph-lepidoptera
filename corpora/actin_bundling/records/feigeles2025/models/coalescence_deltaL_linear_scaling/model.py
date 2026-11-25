def run(inputs, settings):
    # gamma [uN/m] / eta [Pa*s] = (gamma * 1e-6 N/m) / (eta N*s/m^2) = (gamma/eta) um/s
    return {"Delta_L": (inputs["gamma"] / inputs["eta"]) * inputs["tau"]}


def symbolic():
    import sympy as sp

    eta, gamma, tau = sp.symbols("eta gamma tau", positive=True)
    return {"Delta_L": (gamma / eta) * tau}
