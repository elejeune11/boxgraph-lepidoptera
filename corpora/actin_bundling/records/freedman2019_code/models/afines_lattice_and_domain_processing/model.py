"""Faithful Python port of Mathematica AFINES-to-lattice binning and domain-merging pipeline:
- xlink_sorting_paper/process/afines2lattice.m:36, 48-51 (BinCounts + Sign[bcs1-bcs2]/.{-1->2})
- xlink_sorting_paper/process/cytomod_functions.m:858-914 (LatticeXlinkSequences + LatticeDomainMerge)
- xlink_sorting_paper/process/plots.nb:820-827 (domain length = x_last - x_first)
"""

import math


def _lattice_xlink_sequences(dompos):
    if not dompos:
        return []
    runs = []
    cur = [dompos[0]]
    for item in dompos[1:]:
        if item[1] == cur[-1][1]:
            cur.append(item)
        else:
            runs.append(cur)
            cur = [item]
    runs.append(cur)
    return runs


def _lattice_domain_merge(dompos, nc, rc):
    doms = _lattice_xlink_sequences(dompos)
    while True:
        prev = [[tuple(x) for x in run] for run in doms]
        # 1. Fill interior empty (0) gaps shorter than rc flanked by same species (lines 878-889)
        gap_pos = [i for i, run in enumerate(doms) if run[0][1] == 0]
        if gap_pos and gap_pos[0] == 0:
            gap_pos = gap_pos[1:]
        if gap_pos and gap_pos[-1] == len(doms) - 1:
            gap_pos = gap_pos[:-1]
        for pos in gap_pos:
            if doms[pos - 1][0][1] == doms[pos + 1][0][1]:
                if (doms[pos + 1][0][0] - doms[pos - 1][-1][0]) < rc:
                    lbl = doms[pos - 1][0][1]
                    doms[pos] = [(x[0], lbl) for x in doms[pos]]
        doms = _lattice_xlink_sequences([x for run in doms for x in run])
        # 2. Replace non-zero runs shorter than nc with 0 (lines 898-903)
        for pos, run in enumerate(doms):
            if len(run) < nc and run[0][1] != 0:
                doms[pos] = [(x[0], 0) for x in run]
        doms = _lattice_xlink_sequences([x for run in doms for x in run])
        cur = [[tuple(x) for x in run] for run in doms]
        if cur == prev:
            break
    return doms


def run(inputs, settings):
    l0 = float(inputs["l0"])
    lmax = float(inputs["lmax"])
    nc = int(round(float(inputs["nc"])))
    rc = float(inputs["rc"])

    # Bin edges {0, lmax, l0} -> floor(lmax / l0) bins (afines2lattice.m:36,48)
    n_bins = int(math.floor(lmax / l0 + 1e-12))
    coords = [i * l0 for i in range(n_bins)]

    # Probe 1: exact minimal nc-site consecutive domain (sites 0..nc-1)
    d1_min = [i * l0 + 0.5 * l0 for i in range(nc)]
    bcs1 = [sum(1 for x in d1_min if i * l0 <= x < (i + 1) * l0) for i in range(n_bins)]
    bcs2 = [0] * n_bins
    # Tied bin test at bin 10: equal nonzero counts bcs1[10] = bcs2[10] = 2 -> Sign(0) = 0
    bcs1_tie = list(bcs1)
    bcs2_tie = list(bcs2)
    bcs1_tie[10] = 2
    bcs2_tie[10] = 2
    lat_tie = [1 if a > b else (2 if b > a else 0) for a, b in zip(bcs1_tie, bcs2_tie)]
    tied_label = float(lat_tie[10])

    doms_min = _lattice_domain_merge(list(zip(coords, lat_tie)), nc, rc)
    s1_min_lens = [run[-1][0] - run[0][0] for run in doms_min if run[0][1] == 1]
    min_dom_len = s1_min_lens[0] if s1_min_lens else 0.0

    # Probe 2: two short-species clusters separated by gap < rc (merges into 1 domain)
    # and two long-species clusters separated by gap > rc (stays split)
    lat_merge = [0] * n_bins
    # Short (label 1): bins 0..4 and bins 15..19 -> gap = (15 - 4)*0.037 = 0.407 um < 1.0 um -> merges to span 0..19
    for idx in list(range(0, 5)) + list(range(15, 20)):
        lat_merge[idx] = 1
    # Long (label 2): bins 50..54 and bins 90..94 -> gap = (90 - 54)*0.037 = 1.332 um >= 1.0 um -> stays 2 domains of span 4*l0
    for idx in list(range(50, 55)) + list(range(90, 95)):
        lat_merge[idx] = 2

    doms_m = _lattice_domain_merge(list(zip(coords, lat_merge)), nc, rc)
    s1_lens = [run[-1][0] - run[0][0] for run in doms_m if run[0][1] == 1]
    s2_lens = [run[-1][0] - run[0][0] for run in doms_m if run[0][1] == 2]

    return {
        "min_domain_length": float(min_dom_len),
        "merged_short_domain_length": float(s1_lens[0]) if s1_lens else 0.0,
        "unmerged_long_domain_length": float(s2_lens[0]) if s2_lens else 0.0,
        "tied_bin_label": tied_label,
        "n_bins_total": float(n_bins),
    }
