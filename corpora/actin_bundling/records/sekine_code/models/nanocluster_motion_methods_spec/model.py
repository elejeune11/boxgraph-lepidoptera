"""Implementation of the documented nanocluster tracking algorithm described in prose in
README.md (lines 5-8) and sekine2024 Methods p. 10 ('Motion anisotropy'), for direct comparison
against the ported MATLAB implementation (nanocluster_motion_matlab):
  1. Otsu binarization (Methods p. 10; handles bimodal grayscale images with dim clusters).
  2. Filter t_n clusters by physical area A_px * dx^2 in [0.01, 0.1] um^2 AND circularity
     C = 4*pi*A/P^2 in [0.5, 1.0]. No undisclosed border clearing.
  3. Match directly by pixel mask overlap (>= 1 pixel) AND area ratio A2/A1 in [0.5, 2.0]
     (filtering candidate area ratios BEFORE breaking multi-overlap ties by largest overlap,
     and avoiding Stage 1 centroid-miss and Stage 2 hole-shrink bugs).
  4. Convert centroid displacements to nm using pixel_size and compute <|dx|>, <|dy|>,
     mean Euclidean step distance, and anisotropy ratio <|dy|>/<|dx|>.
"""

import math
import numpy as np
from scipy import ndimage

CONN4 = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=np.int32)
CONN8 = np.ones((3, 3), dtype=np.int32)


def otsu_binarize(img):
    """Exact 2-class Otsu thresholding for uint8/float images (reduces to >0 for binary images)."""
    arr = np.asarray(img, dtype=np.float64)
    if arr.min() == arr.max():
        return np.zeros(arr.shape, dtype=bool)
    vals = np.unique(arr)
    if vals.size == 2:
        return arr > vals[0]
    hist, edges = np.histogram(arr.ravel(), bins=256)
    centers = 0.5 * (edges[:-1] + edges[1:])
    w1 = np.cumsum(hist)
    w2 = np.cumsum(hist[::-1])[::-1]
    mean1 = np.cumsum(hist * centers) / np.maximum(w1, 1e-12)
    mean2 = (np.cumsum((hist * centers)[::-1]) / np.maximum(w2[::-1], 1e-12))[::-1]
    var_between = w1[:-1] * w2[1:] * (mean1[:-1] - mean2[1:]) ** 2
    idx = int(np.argmax(var_between))
    return arr > centers[idx]


def edge_perimeter(mask):
    m = np.asarray(mask, dtype=np.int32)
    p = np.pad(m, 1, mode="constant", constant_values=0)
    return float(np.sum(np.abs(np.diff(p, axis=0))) + np.sum(np.abs(np.diff(p, axis=1))))


def circularity_from_mask(mask):
    area = float(np.asarray(mask, dtype=bool).sum())
    perim = edge_perimeter(mask)
    if area <= 0.0 or perim <= 0.0:
        return 0.0
    return min(1.0, 4.0 * math.pi * area / (perim * perim))


def build_scenario(name):
    """Identical deterministic 50x50 synthetic scenarios as nanocluster_motion_matlab."""
    h, w = 50, 50
    img1 = np.zeros((h, w), dtype=np.uint8)
    img2 = np.zeros((h, w), dtype=np.uint8)

    if name == "anisotropic_translation":
        img1[10:14, 10:14] = 255
        img2[12:16, 11:15] = 255
        img1[30:34, 30:34] = 255
        img2[32:36, 29:33] = 255
    elif name == "elongated_low_circularity":
        img1[20:21, 10:30] = 255
        img2[20:21, 12:32] = 255
    elif name == "area_ratio_violation":
        img1[20:24, 20:24] = 255
        img2[20:28, 20:28] = 255
    elif name == "multi_overlap_largest_vs_valid":
        img1[20:24, 20:24] = 255
        img2[22:26, 20:24] = 255
        img2[13:21, 23:31] = 255
    elif name == "border_cluster_removal":
        img1[0:4, 10:14] = 255
        img2[2:6, 11:15] = 255
        img1[44:48, 30:34] = 255
        img2[46:50, 31:35] = 255
    elif name == "antidiagonal_centroid_miss":
        img1[10:13, 10:13] = 255
        img1[12, 12] = 0
        img2[11:15, 11:15] = 255
        img2[11, 11] = 0
        img2[12, 12] = 0
    elif name == "ring_hole_try_catch_abort":
        img1[20:24, 4:8] = 255
        img2[22:26, 5:9] = 255
        img1[20:23, 14:17] = 255
        img1[20:23, 21:23] = 255
        img2[20:23, 14:17] = 255
        img2[21, 15] = 0
        img2[21, 17:21] = 255
        img2[20:23, 21:23] = 255
        img1[20:24, 34:38] = 255
        img2[22:26, 35:39] = 255
    elif name == "grayscale_im2bw_threshold":
        img1[10:14, 10:14] = 100
        img2[12:16, 11:15] = 100
        img1[30:34, 30:34] = 100
        img2[32:36, 29:33] = 100
    else:
        raise ValueError(f"unknown scenario: {name!r}")

    return img1, img2


def track_frame_pair_spec(
    image_1,
    image_2,
    area_min_um2=0.01,
    area_max_um2=0.1,
    circ_min=0.5,
    circ_max=1.0,
    ratio_min=0.5,
    ratio_max=2.0,
    overlap_min_px=1,
    conn=8,
    pixel_size_um=0.04325,
):
    struct = CONN8 if int(conn) == 8 else CONN4
    bw1 = otsu_binarize(image_1)
    bw2 = otsu_binarize(image_2)
    lbl1, n1 = ndimage.label(bw1, structure=struct)
    lbl2, n2 = ndimage.label(bw2, structure=struct)

    a_px_um2 = float(pixel_size_um) ** 2
    px_nm = float(pixel_size_um) * 1000.0

    areas2 = np.bincount(lbl2.ravel()) if n2 > 0 else np.zeros(1, dtype=np.int64)
    centroids2 = {}
    for k in range(1, n2 + 1):
        r2, c2 = np.nonzero(lbl2 == k)
        centroids2[k] = (float(np.mean(c2) + 1.0), float(np.mean(r2) + 1.0))

    matched = []
    unmatched = 0
    dropped = 0

    for k1 in range(1, n1 + 1):
        m1 = lbl1 == k1
        a1 = int(m1.sum())
        a1_um2 = a1 * a_px_um2
        if not (area_min_um2 <= a1_um2 <= area_max_um2):
            continue
        c1_val = circularity_from_mask(m1)
        if not (circ_min <= c1_val <= circ_max):
            dropped += 1
            continue

        r1, c1 = np.nonzero(m1)
        x1, y1 = float(np.mean(c1) + 1.0), float(np.mean(r1) + 1.0)

        # Exact pixel overlap with all t_{n+1} objects
        ov_labels = lbl2[m1]
        ov_labels = ov_labels[ov_labels > 0]
        if ov_labels.size == 0:
            unmatched += 1
            continue

        uniq, counts = np.unique(ov_labels, return_counts=True)
        # Filter candidates by overlap >= overlap_min_px AND area ratio in [ratio_min, ratio_max]
        valid_cands = []
        for k2, ov in zip(uniq.tolist(), counts.tolist()):
            if ov < overlap_min_px:
                continue
            a2 = int(areas2[k2])
            ratio = a2 / float(a1)
            if ratio_min <= ratio <= ratio_max:
                x2, y2 = centroids2[k2]
                dist_px = math.hypot(x2 - x1, y2 - y1)
                valid_cands.append((-ov, dist_px, k2, x2, y2, ratio))

        if not valid_cands:
            unmatched += 1
            continue

        valid_cands.sort()
        _, _, _, x2, y2, ratio = valid_cands[0]
        dx_nm = abs(x2 - x1) * px_nm
        dy_nm = abs(y2 - y1) * px_nm
        matched.append((dx_nm, dy_nm, math.hypot(dx_nm, dy_nm), ratio))

    if not matched:
        return {
            "matched_count": 0.0,
            "unmatched_count": float(unmatched),
            "dropped_count": float(dropped),
            "mean_abs_dx": 0.0,
            "mean_abs_dy": 0.0,
            "mean_step_distance": 0.0,
            "anisotropy_ratio": 0.0,
            "mean_area_ratio": 0.0,
        }

    mean_dx = sum(m[0] for m in matched) / len(matched)
    mean_dy = sum(m[1] for m in matched) / len(matched)
    mean_dist = sum(m[2] for m in matched) / len(matched)
    mean_ratio = sum(m[3] for m in matched) / len(matched)
    aniso = (mean_dy / mean_dx) if mean_dx > 0.0 else 0.0
    return {
        "matched_count": float(len(matched)),
        "unmatched_count": float(unmatched),
        "dropped_count": float(dropped),
        "mean_abs_dx": float(mean_dx),
        "mean_abs_dy": float(mean_dy),
        "mean_step_distance": float(mean_dist),
        "anisotropy_ratio": float(aniso),
        "mean_area_ratio": float(mean_ratio),
    }


def run(inputs, settings):
    scenario = settings.get("scenario", "anisotropic_translation")
    img1, img2 = build_scenario(scenario)
    return track_frame_pair_spec(
        img1,
        img2,
        area_min_um2=float(inputs["area_min"]),
        area_max_um2=float(inputs["area_max"]),
        circ_min=float(inputs["circ_min"]),
        circ_max=float(inputs["circ_max"]),
        ratio_min=float(inputs["ratio_min"]),
        ratio_max=float(inputs["ratio_max"]),
        overlap_min_px=int(round(float(inputs["overlap_min_px"]))),
        conn=int(round(float(inputs["connectivity"]))),
        pixel_size_um=float(inputs["pixel_size"]),
    )
