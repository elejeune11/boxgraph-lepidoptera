"""Shared calibrated reader for Fig. S3 (condensate spreading aspect ratio vs time for Ri = 1.4, 2.2, 2.4 um)."""

from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
from scipy import ndimage
from scipy.signal import fftconvolve

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _calib as C


def calibrate(image: np.ndarray):
    dark, gray, hl, vl = C.find_box_spines(image, thresh=140.0)
    ht, hb = hl[0], hl[-1]
    vl0, vr0 = vl[0], vl[-1]
    # Top 1.1 tick sits 8 px below the top spine, so search starting at ht[2] + 1
    yt_r = C.match_ticks_regular(
        C.find_ticks(dark, vr0, (int(ht[2]) + 1, int(hb[1]) - 1), "y", side=-1, band=(8, 35)), 5
    )
    # Left spine is at ~-5 s; calibrate x on ticks 0..60 plus right spine 70
    xt = C.match_ticks_regular(
        C.find_ticks(dark, hb, (int(vl0[2]) + 15, int(vr0[1]) - 15), "x", side=-1, band=(8, 35)), 7
    )
    xa = C.verify_axis(
        "linear", list(zip(xt, [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0])) + [(vr0[0], 70.0)], 3.0
    )
    ya = C.verify_axis(
        "linear", list(zip(yt_r, [1.1, 1.0, 0.9, 0.8, 0.7])) + [(hb[0], 0.6)], 3.0
    )
    return xa, ya, ht, hb, vl0, vr0


def _make_isolated_template(mask: np.ndarray) -> np.ndarray:
    lab, _ = ndimage.label(mask)
    cands = []
    for i, sl in enumerate(ndimage.find_objects(lab), 1):
        if sl is None:
            continue
        sub = lab[sl] == i
        h, w = sub.shape
        a = int(sub.sum())
        if 35 <= h <= 75 and 35 <= w <= 75 and abs(h - w) <= 8 and a >= 400:
            cands.append((a, sub))
    med_a = float(np.median([c[0] for c in cands]))
    return min(cands, key=lambda c: abs(c[0] - med_a))[1]


def _split_with_template(
    mask: np.ndarray,
    tpl: np.ndarray,
    min_new_frac: float = 0.20,
    min_sep_px: float = 15.0,
    min_xor_gain_frac: float = 0.05,
    max_per_cluster: int = 15,
) -> list[tuple[float, float, int]]:
    tpl_b = tpl.astype(bool)
    tpl_f = tpl_b.astype(np.float32)
    single_area = float(tpl_b.sum())
    th, tw = tpl_b.shape
    cy_t, cx_t = th // 2, tw // 2

    def paint_union(shape: tuple[int, int], pts: list[list[float]]) -> np.ndarray:
        H, W = shape
        canvas = np.zeros((H, W), bool)
        for cr, cc in pts:
            r0 = int(round(cr)) - cy_t
            c0 = int(round(cc)) - cx_t
            sr0, sr1 = max(0, r0), min(H, r0 + th)
            sc0, sc1 = max(0, c0), min(W, c0 + tw)
            if sr0 < sr1 and sc0 < sc1:
                canvas[sr0:sr1, sc0:sc1] |= tpl_b[sr0 - r0 : sr1 - r0, sc0 - c0 : sc1 - c0]
        return canvas

    def refine(padded: np.ndarray, placed: list[list[float]]) -> int:
        H, W = padded.shape
        cov_cnt = np.zeros((H, W), dtype=np.int16)

        def add_tpl(buf: np.ndarray, cr: float, cc: float, delta: int):
            r0 = int(round(cr)) - cy_t
            c0 = int(round(cc)) - cx_t
            sr0, sr1 = max(0, r0), min(H, r0 + th)
            sc0, sc1 = max(0, c0), min(W, c0 + tw)
            if sr0 < sr1 and sc0 < sc1:
                buf[sr0:sr1, sc0:sc1] += delta * tpl_b[sr0 - r0 : sr1 - r0, sc0 - c0 : sc1 - c0]

        for cr, cc in placed:
            add_tpl(cov_cnt, cr, cc, +1)

        for step in (2, 1):
            for _it in range(8):
                improved = False
                for k in range(len(placed)):
                    cr, cc = placed[k]
                    add_tpl(cov_cnt, cr, cc, -1)
                    br0 = max(0, int(round(cr)) - cy_t - step - 1)
                    br1 = min(H, int(round(cr)) - cy_t + th + step + 2)
                    bc0 = max(0, int(round(cc)) - cx_t - step - 1)
                    bc1 = min(W, int(round(cc)) - cx_t + tw + step + 2)
                    other = cov_cnt[br0:br1, bc0:bc1] > 0
                    sub_p = padded[br0:br1, bc0:bc1]

                    def eval_at(tr: float, tc: float) -> int:
                        r0 = int(round(tr)) - cy_t - br0
                        c0 = int(round(tc)) - cx_t - bc0
                        sr0, sr1 = max(0, r0), min(br1 - br0, r0 + th)
                        sc0, sc1 = max(0, c0), min(bc1 - bc0, c0 + tw)
                        comb = other.copy()
                        if sr0 < sr1 and sc0 < sc1:
                            comb[sr0:sr1, sc0:sc1] |= tpl_b[sr0 - r0 : sr1 - r0, sc0 - c0 : sc1 - c0]
                        return int((comb ^ sub_p).sum())

                    best_local = eval_at(cr, cc)
                    best_pos = (cr, cc)
                    for dr, dc in (
                        (step, 0), (-step, 0), (0, step), (0, -step),
                        (step, step), (step, -step), (-step, step), (-step, -step),
                    ):
                        err_l = eval_at(cr + dr, cc + dc)
                        if err_l < best_local:
                            best_local = err_l
                            best_pos = (cr + dr, cc + dc)
                    add_tpl(cov_cnt, best_pos[0], best_pos[1], +1)
                    if best_pos != (cr, cc):
                        placed[k] = [best_pos[0], best_pos[1]]
                        improved = True
                if not improved:
                    break
        return int(((cov_cnt > 0) ^ padded).sum())

    lab, _ = ndimage.label(mask)
    centers: list[tuple[float, float, int]] = []
    for i, sl in enumerate(ndimage.find_objects(lab), 1):
        if sl is None:
            continue
        sub = lab[sl] == i
        area = int(sub.sum())
        if area < 0.25 * single_area:
            continue
        r0, c0 = sl[0].start, sl[1].start
        h, w = sub.shape
        if area <= 1.32 * single_area and h <= 1.25 * th and w <= 1.25 * tw:
            rr, cc = np.nonzero(sub)
            centers.append((float(cc.mean()) + c0 + 0.5, float(rr.mean()) + r0 + 0.5, 1))
            continue

        pad = max(th, tw) + 4
        padded = np.pad(sub, pad, mode="constant", constant_values=False)
        H, W = padded.shape
        Y, X = np.ogrid[:H, :W]
        outside = (~padded).astype(np.float32)
        spill_map = fftconvolve(outside, tpl_f[::-1, ::-1], mode="same")
        placed: list[list[float]] = []
        covered = np.zeros_like(padded)
        cur_err = int(padded.sum())

        for _ in range(max_per_cluster):
            uncov = (padded & ~covered).astype(np.float32)
            if uncov.sum() < 0.15 * single_area:
                break
            gain = fftconvolve(uncov, tpl_f[::-1, ::-1], mode="same")
            score = gain - 1.3 * spill_map
            for pr, pc in placed:
                score[((Y - pr) ** 2 + (X - pc) ** 2) < min_sep_px**2] = -1e9
            br, bc = divmod(int(np.argmax(score)), W)
            if gain[br, bc] < min_new_frac * single_area or spill_map[br, bc] > 0.45 * single_area:
                break
            cand = [list(p) for p in placed] + [[float(br), float(bc)]]
            new_err = refine(padded, cand)
            if cur_err - new_err < min_xor_gain_frac * single_area:
                break
            placed, cur_err = cand, new_err
            covered = paint_union((H, W), placed)

        for pr, pc in placed:
            centers.append((pc - pad + c0 + 0.5, pr - pad + r0 + 0.5, len(placed)))
    return centers


def read_series(image: np.ndarray, which: str) -> list[dict]:
    xa, ya, ht, hb, vl0, vr0 = calibrate(image)
    R = image[:, :, 0].astype(float)
    G = image[:, :, 1].astype(float)
    B = image[:, :, 2].astype(float)
    inb = np.zeros(image.shape[:2], bool)
    inb[int(ht[2]) + 25 : int(hb[1]) - 25, int(vl0[2]) + 25 : int(vr0[1]) - 25] = True
    # Exclude legend box (x > 42, y > 0.95) and droplet cartoon (x in [20, 44], y > 0.97)
    inb[: int(ya.to_pixel(0.95)), int(xa.to_pixel(42.0)) :] = False
    inb[: int(ya.to_pixel(0.97)), int(xa.to_pixel(20.0)) : int(xa.to_pixel(44.0))] = False

    if which == "Ri_1p4":
        mask = inb & (R > 130) & (R - G > 50) & (R - B > 50)
        label = "FigS3 red X (Ri=1.4 um)"
    elif which == "Ri_2p2":
        mask = inb & (R < 60) & (G < 60) & (B < 60)
        label = "FigS3 black circle (Ri=2.2 um)"
    elif which == "Ri_2p4":
        mask = inb & (B > 100) & (B - R > 40) & (G < 130)
        label = "FigS3 blue diamond (Ri=2.4 um)"
    else:
        raise ValueError(which)

    tpl = _make_isolated_template(mask)
    centers = _split_with_template(mask, tpl)
    out = []
    for c, r, mult in centers:
        half_px = (3.0 if mult == 1 else 5.5) + ya.residual_px
        out.append(
            C.make_marker_reading(
                c,
                r,
                xa,
                ya,
                half_px,
                f"{label}, cluster={mult}, cal_res=({xa.residual_px:.2f},{ya.residual_px:.2f})px",
            )
        )
    out.sort(key=lambda d: d["x"])
    return out
