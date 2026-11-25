import math


def run(inputs, settings):
    A = inputs["A"]
    B = inputs["B"]
    alpha = inputs["alpha"]
    t = inputs["t"]
    a = (A * t + B) ** alpha
    a_0 = B ** alpha
    t_0 = B / A
    t_plot = t + t_0
    alpha_recovered = math.log(a / a_0) / math.log(t_plot / t_0)
    return {
        "a": float(a),
        "a_0": float(a_0),
        "t_plot": float(t_plot),
        "alpha_recovered": float(alpha_recovered),
    }


def symbolic():
    import sympy as sp
    A, B, alpha, t = sp.symbols("A B alpha t", positive=True)
    a = (A * t + B) ** alpha
    a_0 = B ** alpha
    t_0 = B / A
    t_plot = t + t_0
    alpha_recovered = sp.log(a / a_0) / sp.log(t_plot / t_0)
    return {"a": a, "a_0": a_0, "t_plot": t_plot, "alpha_recovered": alpha_recovered}
