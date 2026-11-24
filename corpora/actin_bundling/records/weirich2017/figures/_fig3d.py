"""Shared calibrated reader for Fig. 3D (relaxation time tau vs capillary velocity gamma_I/eta)."""

from __future__ import annotations

import numpy as np
from scipy import ndimage, signal
import _calib as C


def calibrate(image: np.ndarray):
    dark, gray, hl, vl = C.find_box_spines(image, row_frac=0.5, col_frac=0.5)
    ht, hb = hl[0], hl[-1]
    vl0, vr0 = vl[0], vl[-1]
    xt = C.find_ticks(dark, hb, (int(vl0[2]) + 5, int(vr0[1]) - 5), "x", side=-1, band=(8, 45))
    xt = C.match_ticks_regular(xt, 18, tol=4.0)
    yt = C.find_ticks(dark, vl0, (int(ht[2]) + 5, int(hb[1]) - 5), "y", side=+1, band=(8, 45))
    yt = C.match_ticks_regular(yt, 12, tol=4.0)

    xa = C.verify_axis("linear", [(vl0[0], 0.0)] + list(zip(xt, [10.0 * i for i in range(1, 19)])), 1.0)
    ya = C.verify_axis("linear", list(zip(yt, [60.0 - 5.0 * i for i in range(12)])) + [(hb[0], 0.0)], 1.0)
    return xa, ya, int(ht[2]) + 65, int(hb[1]) - 65, int(vl0[2]) + 65, int(vr0[1]) - 65


def read_series(image: np.ndarray, at: list[float], which: str) -> list[dict]:
    xa, ya, r0, r1, c0, c1 = calibrate(image)
    r = image[:, :, 0].astype(float)
    g = image[:, :, 1].astype(float)
    b = image[:, :, 2].astype(float)

    inb = np.zeros(image.shape[:2], bool)
    inb[r0:r1, c0:c1] = True
    # Exclude legend box in top-right (x > 140, y > 46)
    inb[: int(ya.to_pixel(46.0)), int(xa.to_pixel(140.0)) :] = False

    if which == "10pct":
        blue = inb & (b - r > 50) & (b > 150)
        # Single blue circle has radius ~66 px; distance_transform >= 38 cleanly splits the touching pair at (83, 44) & (90, 44)
        cores = ndimage.distance_transform_edt(blue) >= 38.0
        out = []
        half_px = 3.0 + ya.residual_px
        for bl in C.find_blobs(cores, min_area=300):
            c_m, r_m = bl["c"], bl["r"]
            yv = float(ya.to_value(r_m))
            ya_, yb_ = float(ya.to_value(r_m - half_px)), float(ya.to_value(r_m + half_px))
            out.append(
                {
                    "x": float(xa.to_value(c_m)),
                    "y": yv,
                    "lo": min(ya_, yb_),
                    "hi": max(ya_, yb_),
                    "px": [c_m, r_m],
                    "px_lo": float(r_m - half_px),
                    "px_hi": float(r_m + half_px),
                    "note": f"Fig3D 10% blue circle, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
                }
            )
        out.sort(key=lambda d: (round(d["x"], 1), d["y"]))
        return out

    if which == "5pct":
        red = inb & (r - g > 80) & (r - b > 80)
        # Ignore black dashed theory curve pixels in non-red background penalty
        blk_line = inb & (r < 80) & (g < 80) & (b < 80)
        non_red_bg = (~red) & (~blk_line)
        # Extract empirical cross template from the isolated cross near (100 nm/s, 14.8 s)
        H = 74
        c_iso, r_iso = int(round(xa.to_pixel(100.0))), int(round(ya.to_pixel(14.8)))
        sub_iso = red[r_iso - H : r_iso + H + 1, c_iso - H : c_iso + H + 1]
        tr, tc = np.nonzero(sub_iso)
        cr_c = r_iso - H + int(round(float(tr.mean())))
        cc_c = c_iso - H + int(round(float(tc.mean())))
        tmpl = red[cr_c - H : cr_c + H + 1, cc_c - H : cc_c + H + 1].copy()

        # Greedy matching pursuit at 2x downsample finds all 8 crosses (including the ~23-px stacked pair at 160 nm/s)
        t2 = tmpl[::2, ::2].astype(np.float32)
        r2 = red[::2, ::2].astype(np.float32)
        bg2 = non_red_bg[::2, ::2].astype(np.float32)
        t2_sum = float(t2.sum())
        pen = signal.fftconvolve(bg2, t2[::-1, ::-1], mode="same")
        uncov = r2.copy()
        h2 = t2.shape[0] // 2
        t_dil = ndimage.binary_dilation(t2 > 0.5, iterations=1)
        centers: list[list[int]] = []
        for _ in range(12):
            gain = signal.fftconvolve(uncov, t2[::-1, ::-1], mode="same")
            valid = (pen <= 0.10 * t2_sum) & (gain >= 0.14 * t2_sum)
            if not np.any(valid):
                break
            score = np.where(valid, gain - 2.0 * pen, -1e9)
            ri, ci = np.unravel_index(int(np.argmax(score)), score.shape)
            if score[ri, ci] < 0:
                break
            centers.append([2 * ri, 2 * ci])
            uncov[ri - h2 : ri - h2 + t2.shape[0], ci - h2 : ci - h2 + t2.shape[1]] *= ~t_dil

        # Coordinate descent minimizing exact binary XOR(red, union) ignoring black dashes
        def _make_union(ctrs: list[list[int]], skip: int = -1) -> np.ndarray:
            u = np.zeros_like(red)
            for idx, (rc, cc) in enumerate(ctrs):
                if idx == skip:
                    continue
                u[rc - H : rc + H + 1, cc - H : cc + H + 1] |= tmpl
            return u

        for _ in range(2):
            for idx in range(len(centers)):
                base = _make_union(centers, skip=idx)
                best = (10**12, centers[idx][0], centers[idx][1])
                for step_s, rng_r, rng_c in [(2, 16, 6), (1, 3, 3)]:
                    rc0, cc0 = best[1], best[2]
                    for dr in range(-rng_r, rng_r + 1, step_s):
                        for dc in range(-rng_c, rng_c + 1, step_s):
                            rc, cc = rc0 + dr, cc0 + dc
                            sl = (slice(rc - H, rc + H + 1), slice(cc - H, cc + H + 1))
                            newly = (~base[sl]) & tmpl
                            delta = int(np.sum(newly & non_red_bg[sl])) - int(np.sum(newly & red[sl]))
                            if delta < best[0]:
                                best = (delta, rc, cc)
                centers[idx] = [best[1], best[2]]

        out = []
        half_px = 4.0 + ya.residual_px
        for r_c, c_c in centers:
            c_m, r_m = float(c_c) + 0.5, float(r_c) + 0.5
            yv = float(ya.to_value(r_m))
            ya_, yb_ = float(ya.to_value(r_m - half_px)), float(ya.to_value(r_m + half_px))
            out.append(
                {
                    "x": float(xa.to_value(c_m)),
                    "y": yv,
                    "lo": min(ya_, yb_),
                    "hi": max(ya_, yb_),
                    "px": [c_m, r_m],
                    "px_lo": float(r_m - half_px),
                    "px_hi": float(r_m + half_px),
                    "note": f"Fig3D 5% red cross, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
                }
            )
        out.sort(key=lambda d: (round(d["x"], 1), d["y"]))
        return out

    if which in ("w2p0", "w1p4"):
        # Pure neutral black dashes; inset past inward tick tips (65 px) and mask text labels 'w=2.0' and 'w=1.4'
        blk = inb & (r < 45) & (g < 45) & (b < 45) & (np.abs(r - g) < 15) & (np.abs(r - b) < 15)
        blk[int(ya.to_pixel(25.5)) : int(ya.to_pixel(18.0)), int(xa.to_pixel(68.0)) : int(xa.to_pixel(106.0))] = False
        blk[int(ya.to_pixel(36.0)) : int(ya.to_pixel(29.0)), int(xa.to_pixel(125.0)) : int(xa.to_pixel(165.0))] = False
        # Separate the two hyperbolas cleanly by product (gamma_I/eta)*tau: w=2.0 (~1278 nm) vs w=1.4 (~3607 nm)
        lab, _ = ndimage.label(blk)
        curve_mask = np.zeros_like(blk)
        for bl in C.find_blobs(blk, min_area=120, max_area=3000):
            xv = float(xa.to_value(bl["c"]))
            yv = float(ya.to_value(bl["r"]))
            prod = xv * yv
            if (which == "w2p0" and 900.0 <= prod <= 1800.0) or (which == "w1p4" and 2800.0 <= prod <= 4400.0):
                sl = (slice(bl["r0"], bl["r1"] + 1), slice(bl["c0"], bl["c1"] + 1))
                curve_mask[sl] |= lab[sl] == bl["label"]
        return [
            C.sample_dashed_curve_y(
                curve_mask, xa, ya, float(x), r0, r1, max_search_dx=220, sigma_dx=85.0, deg=2, note_prefix=f"Fig3D {which}: "
            )
            for x in at
        ]

    raise ValueError(which)
