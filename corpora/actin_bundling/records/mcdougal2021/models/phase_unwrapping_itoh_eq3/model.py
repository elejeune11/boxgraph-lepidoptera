"""1D Itoh phase unwrapping algorithm (Eq. 3, p. 7 & SI p. 11; Itoh 1982, ref-46).

phi_unwrapped(m) = phi_wrapped(0) + sum_{n=1}^{m} W(phi_wrapped(n) - phi_wrapped(n - 1)),
where W maps its argument into (-pi, pi].
"""

import math


def _wrap(phi: float) -> float:
    k = math.floor((math.pi - phi) / (2.0 * math.pi))
    w = phi + 2.0 * math.pi * k
    if w <= -math.pi:
        w += 2.0 * math.pi
    elif w > math.pi:
        w -= 2.0 * math.pi
    return w


def run(inputs, settings):
    curve = inputs["phi_wrapped"]
    xs = [float(v) for v in curve["x"]]
    ys = [float(v) for v in curve["y"]]
    if not ys:
        return {"phi_unwrapped": {"x": [], "y": []}}
    unwrapped = [ys[0]]
    acc = ys[0]
    for n in range(1, len(ys)):
        acc += _wrap(ys[n] - ys[n - 1])
        unwrapped.append(acc)
    return {"phi_unwrapped": {"x": xs, "y": unwrapped}}
