import math
from scipy.optimize import minimize_scalar


def run(inputs, settings):
    p1m = float(inputs["phi1_mid"])
    p2m = float(inputs["phi2_mid"])
    c1 = float(inputs["c1"])
    c2 = float(inputs["c2"])
    chi = float(inputs["chi"])

    top_obs = inputs["ci_top_obs"]
    bot_obs = inputs["ci_bot_obs"]
    ru_obs = inputs["ci_right_upper_obs"]
    arm_obs = inputs["blue_arm_obs"]

    def eval_all(al):
        s = math.sqrt(max(0.0, al * al - chi * chi)) / 6.0
        off_v = (c2 * c2 - c1 * c1) / 3.0 + s
        off_h = (c1 * c1 - c2 * c2) / 3.0 + s

        y_top = [p2m + math.sqrt(max(1e-12, (x - p1m) ** 2 + off_v)) for x in top_obs["x"]]
        y_bot = [p2m - math.sqrt(max(1e-12, (x - p1m) ** 2 + off_v)) for x in bot_obs["x"]]
        y_ru = [p2m + math.sqrt(max(1e-12, (x - p1m) ** 2 - off_h)) for x in ru_obs["x"]]
        y_arm = []
        for x in arm_obs["x"]:
            dx = x - p1m
            d11 = 4.0 * c1 * c1 - 12.0 * dx * dx
            d22 = (chi * chi - al * al) / d11
            dy2 = (4.0 * c2 * c2 - d22) / 12.0
            y_arm.append(p2m + math.sqrt(max(1e-12, dy2)))
        return y_top, y_bot, y_ru, y_arm

    def loss(al):
        yt, yb, yru, yarm = eval_all(al)
        err = 0.0
        for pred, obs in ((yt, top_obs["y"]), (yb, bot_obs["y"]), (yru, ru_obs["y"]), (yarm, arm_obs["y"])):
            for p, o in zip(pred, obs):
                err += (p - o) ** 2
        return err

    res = minimize_scalar(loss, bounds=(abs(chi) + 0.01, 0.60), method="bounded")
    al_fit = float(res.x)
    yt, yb, yru, yarm = eval_all(al_fit)
    return {
        "alpha_fit": al_fit,
        "ci_top_fit": {"x": list(top_obs["x"]), "y": yt},
        "ci_bot_fit": {"x": list(bot_obs["x"]), "y": yb},
        "ci_right_upper_fit": {"x": list(ru_obs["x"]), "y": yru},
        "blue_arm_fit": {"x": list(arm_obs["x"]), "y": yarm},
    }
