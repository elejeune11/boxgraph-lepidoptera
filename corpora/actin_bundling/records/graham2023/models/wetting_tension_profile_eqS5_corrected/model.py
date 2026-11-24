import math


def run(inputs, settings):
    gamma_AVR = inputs["gamma_AVR"]
    gamma_AVD = inputs["gamma_AVD"]
    e = inputs["e"]
    k_e = inputs["k_e"]
    Delta = 0.5 * (gamma_AVD - gamma_AVR)
    return {"gamma": Delta * math.tanh(k_e * (e - 1.0)) + Delta + gamma_AVR}


def symbolic():
    import sympy as sp

    gamma_AVR, gamma_AVD, e, k_e = sp.symbols("gamma_AVR gamma_AVD e k_e", positive=True)
    Delta = sp.Rational(1, 2) * (gamma_AVD - gamma_AVR)
    return {"gamma": Delta * sp.tanh(k_e * (e - 1)) + Delta + gamma_AVR}
