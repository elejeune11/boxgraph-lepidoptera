import math
import numpy as np


def run(inputs, settings):
    sf = (inputs["omega"] / inputs["omega_ref"]) ** 2
    Rav = inputs["C_actin"] / inputs["C_vasp"]
    sigma_o = inputs["sigma_o"]
    sigma_a = inputs["sigma_a"]
    Ka = inputs["K_a"]
    beta = inputs["beta"]
    Xm = inputs["X_min"]
    Xb = inputs["X_max"]

    s = [sf * (sigma_o + i * sigma_a) for i in range(5)]
    sa = sf * sigma_a

    def integral(sig):
        if abs(sig) < 1e-12:
            return Xb - Xm
        return (math.exp(sig * Xb) - math.exp(sig * Xm)) / sig

    I = [integral(si) for si in s]
    Ia = integral(sa)

    B0 = -Rav * I[0]
    B1 = (1.0 - Rav) * Ka * I[1] + Ia * I[0]
    B2 = (2.0 - Rav) * beta * (Ka ** 2) * I[2] + Ka * Ia * I[1]
    B3 = (3.0 - Rav) * (beta ** 3) * (Ka ** 3) * I[3] + beta * (Ka ** 2) * Ia * I[2]
    B4 = (4.0 - Rav) * (beta ** 6) * (Ka ** 4) * I[4] + (beta ** 3) * (Ka ** 3) * Ia * I[3]
    B5 = (beta ** 6) * (Ka ** 4) * Ia * I[4]

    roots = np.roots([B5, B4, B3, B2, B1, B0])
    pos_real = [float(r.real) for r in roots if abs(r.imag) < 1e-9 and r.real > 0.0]
    if not pos_real:
        raise ValueError("No positive real root found for A_o in Fig. S5 quintic equation")
    Ao = pos_real[0]

    # Eqn. 5 (Fig. S5 p. 22)
    denom = I[0] + Ka * Ao * (
        I[1] + beta * Ka * Ao * (
            I[2] + (beta ** 2) * Ka * Ao * (
                I[3] + (beta ** 3) * Ka * Ao * I[4]
            )
        )
    )
    Vo = 1.0 / denom

    coop_powers = [0, 0, 1, 3, 6]
    coeffs = [Vo * (beta ** coop_powers[i]) * ((Ka * Ao) ** i) for i in range(5)]
    species_fracs = [coeffs[i] * I[i] for i in range(5)]
    mean_bound = sum(i * species_fracs[i] for i in range(5))

    def V_of_X(X):
        return sum(coeffs[i] * math.exp(s[i] * X) for i in range(5))

    return {
        "A_o": Ao,
        "V_o": Vo,
        "mean_actin_bound": mean_bound,
        "fractional_saturation": mean_bound / 4.0,
        "V_base_to_meniscus": V_of_X(Xb) / V_of_X(Xm),
    }
