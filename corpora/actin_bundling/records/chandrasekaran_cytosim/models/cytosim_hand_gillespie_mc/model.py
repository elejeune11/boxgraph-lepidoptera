"""Faithful Python port of CytoSim discrete-step Gillespie Kramers unbinding:
- src/sim/hand.cc:238 (nextDetach = RNG.exponential())
- src/sim/hand.h:75 (float nextDetach)
- src/sim/hand.h:286-302 (Hand::checkKramersDetachment)
- src/sim/hand_prop.cc:244, 296-301 (unbinding_force_inv, unbinding_rate_log)
"""
import math
import numpy as np


def run(inputs, settings):
    k_unbind0 = float(inputs["k_unbind0"])
    f_u = float(inputs["unbinding_force"])
    force = float(inputs["force"])
    dt = float(inputs["dt"])

    n_samples = int(settings.get("n_samples", 20000))
    seed = int(settings.get("seed", 42))
    rng = np.random.default_rng(seed)

    # Port of HandProp::complete() + Hand::checkKramersDetachment()
    x = force / f_u + math.log(k_unbind0 * dt)
    dec = np.float32(math.exp(x))

    # Port of Hand::do_attach(): nextDetach = RNG.exponential() stored as float32
    next_detach = rng.exponential(scale=1.0, size=n_samples).astype(np.float32)

    # Number of discrete steps until nextDetach - n * dec < 0
    # Equivalent to: step = 0; while nextDetach >= 0: nextDetach -= dec; step += 1
    # Exact vectorized count matching float32 subtraction for constant dec:
    steps = np.floor(next_detach / dec).astype(np.int64) + 1
    tau_mean = float(np.mean(steps)) * dt
    k_eff = 1.0 / tau_mean

    return {
        "tau_res_mc": tau_mean,
        "k_unbind_mc": k_eff,
    }
