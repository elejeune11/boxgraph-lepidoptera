"""Calibrated readers for SI Fig. S7(a) state slice and Fig. S7(d) traveling wave transit period."""

from __future__ import annotations

import numpy as np
from scipy.signal import find_peaks
import _calib as C


def read_figS7a_state_r06(image: np.ndarray, at: list[float]) -> list[dict]:
    """Read Fig. S7(a) dynamical state code along the rho_A^0 = 0.6 horizontal slice."""
    dark, gray, ht, hb, vl, vr = C.find_box(image)
    xt, yt = C.find_inward_ticks(dark, ht, hb, vl, vr)
    x_axis = C.verify_axis("linear", list(zip(C.match_ticks_regular(xt, 4), [50.0, 60.0, 70.0, 80.0])))
    y_axis = C.verify_axis("linear", list(zip(C.match_ticks_regular(yt, 3), [0.8, 0.6, 0.4])))

    r06 = float(y_axis.to_pixel(0.6))
    r06_i = int(round(r06))
    H, W = image.shape[:2]

    out = []
    for v0 in at:
        c_f = float(x_axis.to_pixel(v0))
        c_i = int(round(c_f))
        patch = image[max(0, r06_i - 25) : min(H, r06_i + 26), max(0, c_i - 25) : min(W, c_i + 26)].astype(float)
        n_blue = int(np.sum((patch[:, :, 2] > 180) & (patch[:, :, 0] < 80) & (patch[:, :, 1] < 80)))
        n_red = int(np.sum((patch[:, :, 0] > 180) & (patch[:, :, 1] < 80) & (patch[:, :, 2] < 80)))
        if n_red > 200:
            code, label = 2.0, "traveling (red triangle)"
        elif n_blue > 200:
            code, label = 1.0, "static (blue square)"
        else:
            code, label = 0.0, "homogeneous (gray circle)"
        out.append({
            "x": float(v0),
            "y": code,
            "lo": code - 0.1,
            "hi": code + 0.1,
            "px": [c_f, r06],
            "px_lo": r06 - 14.0,
            "px_hi": r06 + 14.0,
            "note": f"Fig. S7(a) rho_A^0=0.6, v_0={v0:g}: {label} (blue_px={n_blue}, red_px={n_red}, cal_res=({x_axis.residual_px:.2f},{y_axis.residual_px:.2f})px)",
        })
    return out


def read_figS7d_period(image: np.ndarray, at: list[float]) -> list[dict]:
    """Read steady traveling box-transit period T from Fig. S7(d) kymograph at rho_A^0=0.6, v_0=70."""
    dark, gray, ht, hb, vl, _ = C.find_box(image)
    v_lines = C.clusters(np.where(dark.mean(axis=0) > 0.45)[0])
    if len(v_lines) < 2:
        raise RuntimeError(f"Could not detect main kymograph right border in Fig. S7(d): {v_lines}")
    vr = v_lines[1]  # main kymograph box right edge (before colorbar)

    _, yt = C.find_inward_ticks(dark, ht, hb, vl, vr)
    yt_mid = min(yt, key=lambda r: abs(r - 0.5 * (ht[0] + hb[0])))
    t_axis = C.verify_axis("linear", [(ht[0], 600.0), (yt_mid, 300.0), (hb[0], 0.0)], max_residual_px=2.5)

    r0, r1 = ht[2] + 4, hb[1] - 4
    ts = np.array([t_axis.to_value(r + 0.5) for r in range(r0, r1)], dtype=float)
    sel = (ts >= 335.0) & (ts <= 598.0)
    ts_sel = ts[sel]

    periods = []
    ref_px = None
    for frac in (0.20, 0.35, 0.50, 0.65, 0.75, 0.85):
        cx = int(round(vl[0] + frac * (vr[0] - vl[0])))
        strip = image[r0:r1, cx - 5 : cx + 6].astype(float)
        prof = strip[:, :, 0].mean(axis=1) - strip[:, :, 2].mean(axis=1)
        pks, _ = find_peaks(prof[sel], distance=160, prominence=40.0)
        p_ts = np.sort(ts_sel[pks])
        if len(p_ts) >= 2:
            diffs = np.diff(p_ts)
            # Keep clean steady single-period intervals (exclude boundary-clipped or transient intervals)
            good = diffs[(diffs > 65.0) & (diffs < 88.0)]
            periods.extend(good.tolist())
            if abs(frac - 0.75) < 1e-6 and len(pks) >= 2:
                # Anchor overlay cross on a steady crest at x = 3*pi/2
                p_rows = np.arange(r0, r1)[sel][pks]
                ref_px = [float(cx), float(p_rows[1])]

    if len(periods) < 3:
        raise RuntimeError(f"Too few clean inter-crest periods detected in Fig. S7(d): {periods}")

    arr = np.array(periods, dtype=float)
    T_med = float(np.median(arr))
    # Honest 95% interval combining inter-slice spread (2*std + half-range) and t-axis calibration uncertainty
    px_per_t = abs(hb[0] - ht[0]) / 600.0
    cal_dt = (2.0 * max(t_axis.residual_px, 1.0) + 4.0) / px_per_t
    half_w = float(max(2.0 * np.std(arr, ddof=1), 0.5 * (arr.max() - arr.min())) + cal_dt)

    if ref_px is None:
        ref_px = [float(0.5 * (vl[0] + vr[0])), float(t_axis.to_pixel(450.0))]

    out = []
    for v0 in at:
        out.append({
            "x": float(v0),
            "y": T_med,
            "lo": T_med - half_w,
            "hi": T_med + half_w,
            "px": ref_px,
            "px_lo": ref_px[1] - 12.0,
            "px_hi": ref_px[1] + 12.0,
            "note": (
                f"Fig. S7(d) steady crest-to-crest period across {len(periods)} intervals in t in [335,598]: "
                f"median={T_med:.2f}, range=[{arr.min():.2f},{arr.max():.2f}], cal_res={t_axis.residual_px:.2f}px"
            ),
        })
    return out
