"""Faithful line-by-line Python port of Nanocluster_motion.m (lines 28-250).

Source: corpus/actin_bundling/code/sekine_code/Nanocluster_motion.m (Mustafa M. Sami, 2024;
released with Sekine et al., Nature Communications 15, 464, 2024).

Faithfully reproduces every MATLAB Image Processing Toolbox operation and behavioral quirk:
  - Lines 35-36:  im2bw(I) with implicit default normalized threshold level = 0.5
  - Line 40:      bwareafilt(BW_Image_1, [6 53]) (8-connected, inclusive pixel bounds)
  - Lines 41, 45: imclearborder (8-connected removal of ROI boundary-touching clusters in BOTH frames)
  - Line 44:      bwareafilt(BW_Image_2, [1 Inf]) (no upper/ratio area filter in t_{n+1})
  - Lines 50-72:  Stage 1 intersection (immultiply) -> bwlabel(8) -> Centroid (WITHOUT shrink!) -> bwselect(8)
  - Lines 78-242: Stage 2 per-cluster loop inside bare try/catch:
                    * >= 2 intersection components (lines 101-171): imfill('holes') -> bwmorph('shrink', Inf)
                      -> Centroid -> bwselect(Cumm_Img) -> union -> pick LARGEST 8-component (column-major tie break)
                    * == 1 intersection component (lines 173-209): imfill('holes') -> bwmorph('shrink', Inf)
                      -> Centroid -> bwselect(Cumm_Img) -> Area & Centroid
                    * == 0 intersection components (lines 211-236): sentinel (0, 0, 0) in cols 4-6
                    * Any runtime error (e.g., empty Full_Obj from hole-shrink) silently aborts the rest of the frame
  - No circularity filter, no 0.5-2x area-ratio filter, no pixel-to-micron conversion, and no dx/dy/|dy|/|dx|
    calculation exist in Nanocluster_motion.m; compute_displacements_and_anisotropy() provides optional
    downstream physical summary calculations on the exported 6-column table (x1, y1, A1, x2, y2, A2).
"""

import math
import numpy as np
from scipy import ndimage

CONN4 = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=np.int32)
CONN8 = np.ones((3, 3), dtype=np.int32)


def _conn_struct(conn):
    return CONN8 if int(conn) == 8 else CONN4


def _matlab_round(val):
    """MATLAB round(): half away from zero."""
    return int(math.floor(val + 0.5)) if val >= 0.0 else int(math.ceil(val - 0.5))


def matlab_im2bw(img, level=0.5):
    """Port of MATLAB im2bw(I, level) (lines 35-36, 59, 70, 104, 123, 142)."""
    arr = np.asarray(img)
    if arr.dtype == np.bool_:
        return arr.copy()
    if arr.dtype == np.uint8:
        return arr > (level * 255.0)
    if arr.dtype == np.uint16:
        return arr > (level * 65535.0)
    return arr.astype(np.float64) > float(level)


def matlab_bwlabel(bw, conn=8):
    """Port of MATLAB bwlabel(BW, conn) with exact MATLAB column-major component numbering."""
    mask = np.asarray(bw, dtype=bool)
    raw, n = ndimage.label(mask, structure=_conn_struct(conn))
    if n == 0:
        return np.zeros(mask.shape, dtype=np.int32), 0
    # Order labels by first encountered pixel in column-major scan: min (col, row)
    rows, cols = np.nonzero(raw)
    lbls = raw[rows, cols]
    # Pack key = col * H + row so minimum key gives column-major first pixel
    h = mask.shape[0]
    keys = cols.astype(np.int64) * h + rows.astype(np.int64)
    first_key = np.full(n + 1, np.iinfo(np.int64).max, dtype=np.int64)
    np.minimum.at(first_key, lbls, keys)
    order = np.argsort(first_key[1:]) + 1  # old labels sorted by MATLAB column-major order
    remap = np.zeros(n + 1, dtype=np.int32)
    remap[order] = np.arange(1, n + 1, dtype=np.int32)
    return remap[raw], int(n)


def matlab_bwareafilt(bw, area_range, conn=8):
    """Port of MATLAB bwareafilt(BW, [lo hi], conn) (lines 40, 44)."""
    lo, hi = area_range
    labels, n = matlab_bwlabel(bw, conn=conn)
    if n == 0:
        return np.zeros_like(bw, dtype=bool)
    counts = np.bincount(labels.ravel())
    keep = np.zeros(n + 1, dtype=bool)
    for k in range(1, n + 1):
        if lo <= counts[k] <= hi:
            keep[k] = True
    return keep[labels]


def matlab_imclearborder(bw, conn=8):
    """Port of MATLAB imclearborder(BW, conn) (lines 41, 45; default 8-connectivity in 2D)."""
    labels, n = matlab_bwlabel(bw, conn=conn)
    if n == 0:
        return np.zeros_like(bw, dtype=bool)
    border_labels = set(np.unique(labels[0, :]).tolist())
    border_labels.update(np.unique(labels[-1, :]).tolist())
    border_labels.update(np.unique(labels[:, 0]).tolist())
    border_labels.update(np.unique(labels[:, -1]).tolist())
    border_labels.discard(0)
    keep = np.ones(n + 1, dtype=bool)
    keep[0] = False
    for b in border_labels:
        keep[b] = False
    return keep[labels]


def matlab_regionprops_single(bw):
    """Port of MATLAB regionprops(BW, 'Area', 'Centroid') on a single-component logical mask.
    Returns (area, x_1based, y_1based) or raises IndexError (matching MATLAB empty struct access).
    """
    rows, cols = np.nonzero(bw)
    if rows.size == 0:
        raise IndexError("regionprops returned empty struct (Full_Obj is all-zero)")
    # MATLAB Centroid is [x, y] = [mean(col) + 1, mean(row) + 1] (1-based)
    return float(rows.size), float(np.mean(cols) + 1.0), float(np.mean(rows) + 1.0)


def matlab_bwselect(bw, x_1based, y_1based, conn=8):
    """Port of MATLAB bwselect(BW, c, r, n) (lines 67, 119, 188; default n = 8)."""
    mask = np.asarray(bw, dtype=bool)
    h, w = mask.shape
    c0 = _matlab_round(float(x_1based)) - 1
    r0 = _matlab_round(float(y_1based)) - 1
    if not (0 <= r0 < h and 0 <= c0 < w) or not mask[r0, c0]:
        return np.zeros_like(mask, dtype=bool)
    labels, _ = matlab_bwlabel(mask, conn=conn)
    return labels == labels[r0, c0]


def matlab_imfill_holes(bw):
    """Port of MATLAB imfill(BW, 'holes') (lines 110, 181; 4-connected background holes in 2D)."""
    return ndimage.binary_fill_holes(np.asarray(bw, dtype=bool), structure=CONN4)


def _build_shrink_lut():
    """Precompute 512-entry lookup table for 8-connected Euler-preserving shrink (bwmorph 'shrink').
    A center pixel p is removable iff:
      1. p == 1 and has at least one 8-neighbor (isolated pixels stay)
      2. at least one 4-neighbor is 0 (outer boundary peels inward)
      3. p is (8,4)-simple: its 8-neighbors form 1 8-connected foreground component AND its
         4-adjacent background neighbors belong to 1 4-connected background component in the 3x3 patch.
    """
    lut = np.zeros(512, dtype=bool)
    # 8-neighbor adjacency in 3x3 grid excluding center (4)
    fg_coords = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1), (2, 2)]
    fg_adj = {
        i: [j for j in range(8) if i != j and abs(fg_coords[i][0] - fg_coords[j][0]) <= 1
            and abs(fg_coords[i][1] - fg_coords[j][1]) <= 1]
        for i in range(8)
    }
    fg_bit = [0, 1, 2, 3, 5, 6, 7, 8]
    bg_adj4 = {
        i: [j for j in range(8) if i != j
            and (abs(fg_coords[i][0] - fg_coords[j][0]) + abs(fg_coords[i][1] - fg_coords[j][1]) == 1)]
        for i in range(8)
    }
    four_idx = [1, 3, 4, 6]  # positions in fg_coords corresponding to (0,1),(1,0),(1,2),(2,1)

    for code in range(512):
        if not (code & (1 << 4)):
            continue
        # Check at least one 4-neighbor is 0
        if all(code & (1 << fg_bit[i]) for i in four_idx):
            continue
        present = [i for i in range(8) if (code & (1 << fg_bit[i]))]
        if not present:
            continue
        # BFS on 8-connected foreground neighbors
        seen = {present[0]}
        stack = [present[0]]
        while stack:
            u = stack.pop()
            for v in fg_adj[u]:
                if (code & (1 << fg_bit[v])) and v not in seen:
                    seen.add(v)
                    stack.append(v)
        if len(seen) != len(present):
            continue
        # BFS on 4-connected background neighbors among the 8 ring pixels
        bg_four = [i for i in four_idx if not (code & (1 << fg_bit[i]))]
        seen_bg = {bg_four[0]}
        stack = [bg_four[0]]
        while stack:
            u = stack.pop()
            for v in bg_adj4[u]:
                if not (code & (1 << fg_bit[v])) and v not in seen_bg:
                    seen_bg.add(v)
                    stack.append(v)
        if all(i in seen_bg for i in bg_four):
            lut[code] = True
    return lut


_SHRINK_LUT = _build_shrink_lut()
_WEIGHTS_3X3 = np.array([[1, 2, 4], [8, 16, 32], [64, 128, 256]], dtype=np.int32)


def matlab_bwmorph_shrink_inf(bw):
    """Port of MATLAB bwmorph(BW, 'shrink', Inf) (lines 111, 182) via 4 checkerboard subiterations."""
    img = np.asarray(bw, dtype=bool).copy()
    if not img.any():
        return img
    h, w = img.shape
    subgrids = ((0, 0), (1, 1), (0, 1), (1, 0))
    while True:
        changed = False
        for r0, c0 in subgrids:
            padded = np.pad(img, 1, mode="constant", constant_values=False)
            code = np.zeros((h, w), dtype=np.int32)
            for dr in range(3):
                for dc in range(3):
                    code += padded[dr : dr + h, dc : dc + w] * _WEIGHTS_3X3[dr, dc]
            mask = np.zeros((h, w), dtype=bool)
            mask[r0::2, c0::2] = True
            remove = mask & _SHRINK_LUT[code]
            if remove.any():
                img[remove] = False
                changed = True
        if not changed:
            break
    return img


def track_frame_pair_matlab(
    image_1,
    image_2,
    area_min_t1=6,
    area_max_t1=53,
    area_min_t2=1,
    conn=8,
    binarize_level=0.5,
):
    """Faithful port of one iteration of Nanocluster_motion.m (lines 30-247).

    Returns:
        rows: list of [x1, y1, A1, x2, y2, A2] matching ref_Coordinate_List (1-based px coords)
        diag: dict with diagnostic counts (raw_t1_in_size, border_dropped_t1, aborted_by_catch,
              try_catch_dropped_t1)
    """
    # Lines 35-36: im2bw with default threshold 0.5
    bw_1 = matlab_im2bw(image_1, level=binarize_level)
    bw_2 = matlab_im2bw(image_2, level=binarize_level)

    # Lines 40-41: bwareafilt([6 53]) + imclearborder on Image_1
    sized_1 = matlab_bwareafilt(bw_1, (area_min_t1, area_max_t1), conn=conn)
    _, n_sized_1 = matlab_bwlabel(sized_1, conn=conn)
    large_1 = matlab_imclearborder(sized_1, conn=conn)

    # Lines 44-45: bwareafilt([1 Inf]) + imclearborder on Image_2
    sized_2 = matlab_bwareafilt(bw_2, (area_min_t2, float("inf")), conn=conn)
    large_2 = matlab_imclearborder(sized_2, conn=conn)

    # Line 50: immultiply
    intersect_12 = large_1 & large_2

    # Lines 55-72: Stage 1 cumulative overlapped objects in Image_2
    # NOTE: Stage 1 uses Centroid WITHOUT imfill/shrink (non-convex/anti-diagonal centroid-miss bug!)
    labels_int, n_int = matlab_bwlabel(intersect_12, conn=conn)
    cumm_img = np.zeros_like(large_2, dtype=bool)
    for cell_count in range(1, n_int + 1):
        cell_obj = labels_int == cell_count
        _, cx, cy = matlab_regionprops_single(cell_obj)
        full_obj = matlab_bwselect(large_2, cx, cy, conn=conn)
        cumm_img = cumm_img | full_obj

    # Lines 78-84: label ref_Img and Cumm_Img
    ref_labels, ref_n = matlab_bwlabel(large_1, conn=conn)
    border_dropped_t1 = max(0, n_sized_1 - ref_n)

    rows = []
    aborted = False
    # Lines 88-242: bare try ... catch end around the entire per-cluster loop!
    try:
        for ref_cell_count in range(1, ref_n + 1):
            ref_cell_obj = ref_labels == ref_cell_count
            # Line 92: intersection with Cumm_Img
            ref_intersect = ref_cell_obj & cumm_img
            # Lines 96-97: Area of ref_Cell_Obj
            area_ref, ref_x, ref_y = matlab_regionprops_single(ref_cell_obj)
            # Line 99: bwlabel(..., 8)
            labeled_int, n_sub = matlab_bwlabel(ref_intersect, conn=conn)

            if n_sub >= 2:
                # Lines 101-171: >= 2 intersection components -> union matched Cumm_Img objects,
                # then pick LARGEST 8-connected object (maxArea, column-major tie break)
                int_cumm = np.zeros_like(large_2, dtype=bool)
                for many_count in range(1, n_sub + 1):
                    many_obj = labeled_int == many_count
                    many_obj = matlab_imfill_holes(many_obj)
                    many_obj = matlab_bwmorph_shrink_inf(many_obj)
                    _, mx, my = matlab_regionprops_single(many_obj)
                    mmany_obj = matlab_bwselect(cumm_img, mx, my, conn=conn)
                    int_cumm = int_cumm | mmany_obj

                # Lines 130-144: largest overlapped object
                labeled_bw, n_cand = matlab_bwlabel(int_cumm, conn=conn)
                if n_cand == 0:
                    raise IndexError("empty Intersect_Commulative")
                areas = [int(np.sum(labeled_bw == i)) for i in range(1, n_cand + 1)]
                max_area = max(areas)
                max_num = next(i for i, a in enumerate(areas, start=1) if a == max_area)
                full_obj = labeled_bw == max_num
                area_full, fx, fy = matlab_regionprops_single(full_obj)
                rows.append([ref_x, ref_y, area_ref, fx, fy, area_full])

            elif n_sub == 1:
                # Lines 173-209: single intersection component -> imfill + shrink + bwselect
                shrunk = matlab_imfill_holes(ref_intersect)
                shrunk = matlab_bwmorph_shrink_inf(shrunk)
                _, ix, iy = matlab_regionprops_single(shrunk)
                full_obj = matlab_bwselect(cumm_img, ix, iy, conn=conn)
                area_full, fx, fy = matlab_regionprops_single(full_obj)
                rows.append([ref_x, ref_y, area_ref, fx, fy, area_full])

            else:
                # Lines 211-236: no intersection -> put 0 sentinel in columns 4, 5, 6
                rows.append([ref_x, ref_y, area_ref, 0.0, 0.0, 0.0])
    except Exception:
        aborted = True

    try_catch_dropped = max(0, ref_n - len(rows))
    diag = {
        "raw_t1_in_size": n_sized_1,
        "border_dropped_t1": border_dropped_t1,
        "aborted_by_catch": aborted,
        "try_catch_dropped_t1": try_catch_dropped,
    }
    return rows, diag


def build_scenario(name):
    """Deterministic 50x50 ROI synthetic verification scenarios (matching sekine2024 50x50 ROI)."""
    h, w = 50, 50
    img1 = np.zeros((h, w), dtype=np.uint8)
    img2 = np.zeros((h, w), dtype=np.uint8)

    if name == "anisotropic_translation":
        # Two compact 4x4 interior clusters (A=16 px, C=0.785) shifting (+1,+2) and (-1,+2) px
        img1[10:14, 10:14] = 255
        img2[12:16, 11:15] = 255
        img1[30:34, 30:34] = 255
        img2[32:36, 29:33] = 255
    elif name == "elongated_low_circularity":
        # One elongated 1x20 rod (A=20 px in [6,53], C=0.1425 << 0.5) shifting (+2, 0) px
        img1[20:21, 10:30] = 255
        img2[20:21, 12:32] = 255
    elif name == "area_ratio_violation":
        # One 4x4 cluster (A1=16 px) fusing into an 8x8 cable (A2=64 px, ratio=4.0 > 2.0)
        img1[20:24, 20:24] = 255
        img2[20:28, 20:28] = 255
    elif name == "multi_overlap_largest_vs_valid":
        # One 4x4 cluster (A1=16 px) overlapping BOTH a valid 4x4 match (A2=16, ratio=1.0, dy=+2 px)
        # AND an invalid 8x8 cable (A2=64, ratio=4.0, 1-px corner overlap at (20,23))
        img1[20:24, 20:24] = 255
        img2[22:26, 20:24] = 255  # valid match A (rows 22..25, cols 20..23)
        img2[13:21, 23:31] = 255  # invalid large cable B (rows 13..20, cols 23..30; row 21 empty!)
    elif name == "border_cluster_removal":
        # Two compact 4x4 clusters shifting (+1, +2) px:
        # Cluster 1 touches top ROI border (row 0) at t_n; Cluster 2 moves onto bottom border (row 49) at t_{n+1}
        img1[0:4, 10:14] = 255
        img2[2:6, 11:15] = 255
        img1[44:48, 30:34] = 255
        img2[46:50, 31:35] = 255
    elif name == "antidiagonal_centroid_miss":
        # Stage 1 bug (lines 61-67): 2-pixel anti-diagonal overlap {(11,12), (12,11)} whose centroid
        # rounds to (12,12) = background in Image_2 -> Stage 1 bwselect returns empty!
        img1[10:13, 10:13] = 255
        img1[12, 12] = 0  # A1 = 8 px, C = 0.698 > 0.5
        img2[11:15, 11:15] = 255
        img2[11, 11] = 0
        img2[12, 12] = 0  # A2 = 14 px, A2/A1 = 1.75 in [0.5, 2], single 8-connected object
    elif name == "ring_hole_try_catch_abort":
        # Stage 2 bug (lines 88, 181-191, 241-242):
        # Cluster 1 (cols 4..7): healthy 4x4 shifting (+1,+2)
        img1[20:24, 4:8] = 255
        img2[22:26, 5:9] = 255
        # Cluster 2 (cols 14..16): 3x3 solid (A1=9) overlapping a 3x3 hollow ring in t_{n+1} (missing (21,15))
        img1[20:23, 14:17] = 255
        # Cluster 3 (cols 21..22): 3x2 helper (A1=6) overlapping tail of the ring object so Stage 1 adds it to Cumm_Img
        img1[20:23, 21:23] = 255
        # Ring + bridge + tail in t_{n+1} (A2=18 px; ratio 2.0 for Cluster 2, 3.0 for Cluster 3)
        img2[20:23, 14:17] = 255
        img2[21, 15] = 0  # 1-px interior hole!
        img2[21, 17:21] = 255
        img2[20:23, 21:23] = 255
        # Cluster 4 (cols 34..37): completely independent healthy 4x4 shifting (+1,+2) -> SILENTLY DROPPED by try/catch!
        img1[20:24, 34:38] = 255
        img2[22:26, 35:39] = 255
    elif name == "grayscale_im2bw_threshold":
        # Dim uint8 clusters at intensity 100 (100/255 = 0.392 < im2bw default 0.5) over background 0
        img1[10:14, 10:14] = 100
        img2[12:16, 11:15] = 100
        img1[30:34, 30:34] = 100
        img2[32:36, 29:33] = 100
    else:
        raise ValueError(f"unknown scenario: {name!r}")

    return img1, img2


def summarize_rows(rows, diag, pixel_size_um, apply_downstream_area_ratio=False, ratio_min=0.5, ratio_max=2.0):
    """Compute displacement (nm), step distance (nm), and motion anisotropy |dy|/|dx| from ALL_Result rows."""
    px_nm = float(pixel_size_um) * 1000.0
    matched = []
    unmatched = 0
    ratio_rejected = 0

    for x1, y1, a1, x2, y2, a2 in rows:
        if a2 <= 0.0:
            unmatched += 1
            continue
        ratio = a2 / a1
        if apply_downstream_area_ratio and not (ratio_min <= ratio <= ratio_max):
            ratio_rejected += 1
            continue
        dx_nm = abs(x2 - x1) * px_nm
        dy_nm = abs(y2 - y1) * px_nm
        dist_nm = math.hypot(dx_nm, dy_nm)
        matched.append((dx_nm, dy_nm, dist_nm, ratio))

    dropped = int(diag["border_dropped_t1"] + diag["try_catch_dropped_t1"] + ratio_rejected)
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
    apply_ratio = bool(settings.get("apply_downstream_area_ratio", False))
    img1, img2 = build_scenario(scenario)
    rows, diag = track_frame_pair_matlab(
        img1,
        img2,
        area_min_t1=float(inputs["area_min_t1"]),
        area_max_t1=float(inputs["area_max_t1"]),
        area_min_t2=float(inputs["area_min_t2"]),
        conn=int(round(float(inputs["connectivity"]))),
        binarize_level=float(inputs["binarize_level"]),
    )
    return summarize_rows(
        rows,
        diag,
        pixel_size_um=float(inputs["pixel_size"]),
        apply_downstream_area_ratio=apply_ratio,
        ratio_min=float(inputs["ratio_min"]),
        ratio_max=float(inputs["ratio_max"]),
    )
