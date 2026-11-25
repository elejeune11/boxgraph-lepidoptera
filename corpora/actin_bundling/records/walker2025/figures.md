# Figure Reading Log for `walker2025` (Walker et al., *Developmental Cell* 60, 1550–1567, 2025)

All quantitative data plots across Main Figures 1–7 (PDF pp. 4–13, rendered from vector graphics at 1200 dpi) and Supplemental Figures S1–S7 (PDF pp. 29–40, rendered from embedded ~215 ppi rasters at 900 dpi) have been digitized with calibrated `mode: code` readers in `figures.yaml` / `figures/*.py`, linked to datums in `data.yaml`, checked via `boxgraph figure check`, visually verified against overlay PNGs in `figwork/`, and saved to `readings.json`.

---

## Main Figures (PDF pp. 4–13, 1200 dpi vector)

### Figure 1 (p. 4) — *mini-Lpd LLPS*
- **1A–B**: Left out — domain architecture and net-charge-per-residue schematic drawings with no numerical axes/data values.
- **1C–E, 1G**: Left out — fluorescence microscopy images and condensate fusion / FRAP montages with no plotted data.
- **1F**: Read — [`fig1f_frap_minilpd`] (normalized mini-Lpd FRAP recovery ± SD vs time at 13 timestamps from `-6` to `177 s`; `t = 0 s` omitted because it lies on the vertical laser-bleach step). Calibrated from L-axis spines and outward ticks (`x: -50..200 s`, `y: 0.0..1.2`), max calibration residual `x <= 1.2 px`, `y <= 1.4 px`; blue solid curve + black ±SD caps verified on overlay.
- **1H**: Read — [`fig1h_partition`] (mean 10 uM mini-Lpd partition coefficient bar heights across 50, 100, 150, 500 mM NaCl: `49.5, 40.3, 25.7, 2.09`, with overlaid n = 3 replicate circles cross-checked in reading notes). Calibrated from L-axis outward y-ticks (`0..60`), `ya_res <= 1.5 px`.
- **1I**: Read — [`fig1i_diam_50mM`, `fig1i_diam_100mM`, `fig1i_diam_150mM`, `fig1i_diam_500mM`] (droplet diameter frequency polygons across 0–6 um for 50, 100, 150, 500 mM NaCl). Calibrated from L-axis outward ticks (`x: 0..6 um`, `y: 0..500`), `res <= 2.0 px`.

### Figure 2 (p. 5) — *mini-Lpd + VASP co-condensation without PEG*
- **2A, 2I**: Left out — protein domain / interaction schematics with no numerical axes.
- **2B–C, 2E–F, 2H**: Left out — fluorescence micrographs and fusion / FRAP montages.
- **2D**: Read — [`fig2d_diam_3to1`, `fig2d_diam_1to1`, `fig2d_diam_1to3`] (droplet diameter frequency polygons for 3:1, 1:1, 1:3 mini-Lpd:VASP). Calibrated from L-axis outward ticks (`x: 0..6 um`, `y: 0..200`), `res <= 2.0 px`.
- **2G**: Read — [`fig2g_frap_minilpd`, `fig2g_frap_vasp`] (emerald-green mini-Lpd and magenta VASP FRAP recovery % ± SD vs time, n = 9; `t = 0 s` omitted on vertical bleach step). Calibrated from L-axis ticks (`x: -10..190 s`, `y: 0..120%`), `res <= 1.6 px`.

### Figure 3 (p. 6) — *VASP + mini-Lpd actin assembly without PEG*
- **3A–C**: Left out — fluorescence microscopy images and confocal cross-sections.
- **3D**: Left out — schematic mechanism cartoon.
- **3E**: Read — [`fig3e_ar_mean`, `fig3e_ar_median`] (mean black squares `1.025, 1.164, 1.618, 3.008` and horizontal black box medians across `0, 0.25, 0.5, 0.75 uM` actin; IQR box bounds, cyan replicate means, and visible scatter maxima logged in notes). Calibrated from L-axis outward y-ticks (`1..10`), `ya_res <= 1.2 px`.
- **3F**: Read — [`fig3f_frac_ar_gt1p2`] (mean fraction of condensates with AR > 1.2: `0.0107, 0.1479, 0.7076, 0.9247`, n = 3 replicate circles cross-checked). Calibrated from L-axis ticks (`0..1`), `ya_res = 0.00 px`.

### Figure 4 (p. 7) — *mini-Lpd actin assembly & bivalent CytoSim simulations*
- **4A–B, 4E–H**: Left out — fluorescence micrographs and 3D confocal reconstructions.
- **4C**: Read — [`fig4c_ar_mean`, `fig4c_ar_median`, `fig4c_lata_ar_mean`, `fig4c_lata_ar_median`] (mean squares and box medians for `0, 1, 2, 3 uM` actin and `2 uM + LatA` control). Calibrated from L-axis ticks (`1..4`), `ya_res <= 1.5 px`. Note: visible scatter clouds at 2 and 3 uM terminate sharply at `AR = 4.00` without mention in the caption (`q_fig5i_and_ar_cutoffs_verified`).
- **4D**: Read — [`fig4d_frac_ar_gt1p2`, `fig4d_lata_frac_ar_gt1p2`] (mean fraction AR > 1.2: `0.0044, 0.0410, 0.4038, 0.5676` and `0.0036` for `2 uM + LatA`). Calibrated from L-axis ticks (`0..1`), `ya_res = 3.57 px`.
- **4I**: Left out — 3D CytoSim final simulation snapshots (`t = 600 s`) with no plotted axes.
- **4J**: Read — [`fig4j_surf_kb10_ku0p01`, `fig4j_surf_kb10_ku0p1`, `fig4j_surf_kb10_ku1`, `fig4j_surf_kb1_ku0p01`, `fig4j_surf_kb1_ku0p1`, `fig4j_surf_kb1_ku1`, `fig4j_surf_kb0p1_ku0p01`, `fig4j_surf_kb0p1_ku0p1`, `fig4j_surf_kb0p1_ku1`] (all 9 solid mean trajectories of actin-covered surface area fraction across `k_bind = 10, 1.0, 0.1` and `k_unbind = 0.01, 0.1, 1.0 1/s`, `t = 200..590 s`; `t = 300 s` omitted on `kb=0.1, ku=0.1` where grey fully occludes red). Calibrated per stacked box from outward ticks (`x: 0..600 s`, `y: 0.0..0.6`), `res = (0.43, 0.00) px`.
- **4K**: Read — [`fig4k_dimers_kb10_b0..b2`, `fig4k_dimers_kb1_b0..b2`, `fig4k_dimers_kb0p1_b0..b2`] (all 9 stacked-bar bound-filament fraction series 0/1/2 with staggered black ±SD stems across `k_bind = 10, 1.0, 0.1` and `k_unbind = 0.01, 0.1, 1.0 1/s`). Calibrated per box from outward ticks (`y: 0.0, 0.5, 1.0`), `ya_res <= 0.17 px`.

### Figure 5 (p. 9) — *Dynamic dimerization simulations & monomer-mini-Lpd*
- **5A**: Left out — reaction network schematic.
- **5B**: Left out — 3D CytoSim final snapshots (`t = 600 s`).
- **5C**: Read — [`fig5c_surf_kf10_ks0p001`, `fig5c_surf_kf10_ks0p01`, `fig5c_surf_kf10_ks0p1`, `fig5c_surf_kf10_ks1`, `fig5c_surf_kf10_ks10`, `fig5c_surf_kf1_ks0p001`, `fig5c_surf_kf1_ks0p01`, `fig5c_surf_kf1_ks0p1`, `fig5c_surf_kf1_ks10`] (9 cleanly separable solid mean surface-area trajectories across `k_dimer_form = 10, 1.0` and `k_dimer_split = 0.001..10.0 1/s`). Left out curve: `k_dimer_form = 1.0, k_dimer_split = 1.0` (light grey) because it is heavily occluded underneath the `k_dimer_split = 10.0` pale-teal shell trajectory and its SD band. Calibrated per box (`x: 0..600 s`, `y: 0.0..0.6`), `res = (0.54, 0.00) px`.
- **5D**: Read — [`fig5d_kf10_free_monomer`, `fig5d_kf10_bound_monomer`, `fig5d_kf10_free_dimer`, `fig5d_kf10_singly_bound_dimer`, `fig5d_kf10_doubly_bound_dimer`, `fig5d_kf1_free_monomer`, `fig5d_kf1_bound_monomer`, `fig5d_kf1_free_dimer`, `fig5d_kf1_singly_bound_dimer`, `fig5d_kf1_doubly_bound_dimer`] (all 10 stacked-bar population fraction series ± SD across `k_dimer_split = 10.0..0.001 1/s` for both `k_dimer_form = 10` and `1.0 1/s`). Calibrated per box (`y: 0.0..1.0`, 6 ticks), `ya_res <= 0.24 px`. Compared against `well_mixed_dynamic_dimerization_ode` sweeps (`cmp_fig5d_*_ode`, `q_ode_vs_fig5d_spatial_and_cactin`).
- **5E–H**: Left out — fluorescence micrographs and 3D confocal reconstructions.
- **5I**: Read — [`fig5i_ar_mean`, `fig5i_ar_median`] (mean `1.021, 2.507` and median `1.010, 1.112` for monomer-mini-Lpd at `0, 1 uM` actin). Calibrated from L-axis ticks (`1..12`), `ya_res = 2.01 px`. Settles `q_fig5i_axis_vs_caption_cutoff`: visible scatter reaches `max_visible = 12.00` (not clipped at 10 as stated in the caption; `q_fig5i_and_ar_cutoffs_verified`).
- **5J**: Read — [`fig5j_frac_ar_gt1p2`] (mean fraction AR > 1.2: `0.0055, 0.4626`, n = 3). Calibrated from L-axis ticks (`0..1`), `ya_res = 4.07 px`.

### Figure 6 (p. 11) — *WT Eps15 vs 1:1 WT Eps15:Eps15-Lifeact*
- **6A, 6C**: Left out — cartoon schematic + single control micrograph.
- **6B, 6E–I**: Left out — fluorescence micrographs and 3D confocal reconstructions.
- **6D**: Read — [`fig6d_eps15_partition`] (WT Eps15 partition coefficient mean ± SD bars: `12.67` self-partitioning, `12.07` in 1:1 WT:Eps15-Lifeact, n = 3). Calibrated from L-axis ticks (`0..20`), `ya_res = 2.08 px`.
- **6J**: Read — [`fig6j_ar_mean`, `fig6j_ar_median`, `fig6j_wt_ar_mean`, `fig6j_wt_ar_median`] (mean squares `1.015, 1.042, 1.631, 5.003` and medians `1.015, 1.033, 1.250, 4.115` for `0..1.25 uM` actin, plus WT Eps15 control `1.051 / 1.033`). Calibrated from L-axis ticks (`2..16`), `ya_res = 2.27 px`.
- **6K**: Read — [`fig6k_frac_ar_gt1p2`, `fig6k_wt_frac_ar_gt1p2`] (mean fraction AR > 1.2: `0.0060, 0.0304, 0.5784, 0.9620` and WT control `0.0089`). Calibrated on the 6 uniform ticks `0.0..1.0` (`ya_res = 5.65 px`; top tick `1.2` is visually compressed by ~25 px). Note: caption claims ±SD error bars, but none are plotted (`q_fig6k_caption_sd_bars_absent`).

### Figure 7 (p. 13) — *G-actin enrichment & VASP mut-FAB / mut-GAB*
- **7A, 7C, 7E**: Left out — fluorescence microscopy images.
- **7B**: Read — [`fig7b_actin_enrichment`] (broken-y-axis G-actin fold enrichment bar heights: `172.1` mini-Lpd, `27.23` mini-Lpd+VASP, `3.76` Eps15-Lifeact, `3.12` Eps15 WT, directly quantifying the main text's `> 150`, `> 25`, `~3`). Calibrated separately on upper segment (`160..190`, `ya_res = 0.80 px`) and lower segment (`0..30`, `ya_res = 1.50 px`).
- **7D**: Read — [`fig7d_vasp_actin_enrichment`] (G-actin fold enrichment bar heights: `32.31` WT VASP, `22.44` VASP mut-FAB, `19.53` VASP mut-GAB, quantifying text `~20`). Calibrated from L-axis ticks (`0..40`), `ya_res = 4.04 px`.
- **7F**: Read — [`fig7f_ar_mean`, `fig7f_ar_median`] (mean `4.52, 1.069, 1.94` and median `3.285, 1.022, 1.243` for VASP WT, mut-FAB, mut-GAB). Calibrated from L-axis ticks (`2..16`), `ya_res = 1.12 px`. Note: VASP WT scatter terminates sharply at `AR = 15.0` without mention in the caption (`q_fig5i_and_ar_cutoffs_verified`).
- **7G**: Read — [`fig7g_frac_ar_gt1p2`] (mean fraction AR > 1.2: `0.8717, 0.0788, 0.5707`, n = 3). Calibrated from L-axis ticks (`0..1`), `ya_res = 4.08 px`. Note: y-axis misprints `"Fractrion"`.
- **7H**: Left out — schematic mechanism cartoon.

---

## Supplemental Figures (PDF pp. 29–41, 900 dpi raster)

### Figure S1 (pp. 29–30) — *Sequence analysis, monomer partitioning & GUVs*
- **S1A**: Left out — FINCHES 2D residue-pair interaction heatmaps (color-coded 2D matrix, no 1D curve/bar values).
- **S1B**: Read — [`figs1b_minilpd_disorder`, `figs1b_minilpd_plddt`, `figs1b_fllpd_disorder`, `figs1b_fllpd_plddt`] (SPARROW disorder score and predicted pLDDT vs residue index for Lpd aa850–1250 [`25..425`] and full-length Lpd [`50..1200`]). Calibrated from internal light-grey gridlines (`y: 0.0..1.0`, `x: 50..400` and `200..1200`), `res <= (0.96, 2.25) px`.
- **S1C**: Read — [`figs1c_partition_minilpd`, `figs1c_partition_monomer_minilpd`] (mean self-partition coefficient at 50 and 150 mM NaCl: mini-Lpd `49.13, 25.99`, monomer mini-Lpd `40.39, 26.12`, n = 3; mini-Lpd controls closely match Fig. 1H `49.5, 25.7`). Calibrated from L-axis ticks (`0..60`), `ya_res = 3.65 px`.
- **S1D–E**: Left out — fluorescence microscopy images and GUV slices/z-projections.

### Figure S2 (pp. 31–32) — *Condensate deformation time-lapses*
- **S2A–G**: Left out — fluorescence micrographs and sphere-to-rod / elongation microscopy montages with no plotted axes.

### Figure S3 (pp. 33–34) — *FRAP & zippering of actin-deformed condensates*
- **S3A, S3C, S3E, S3G, S3I–K**: Left out — FRAP and condensate-zippering fluorescence microscopy montages.
- **S3B**: Read — [`figs3b_frap_minilpd`, `figs3b_frap_actin`] (normalized FRAP intensity ± SD for mini-Lpd and actin in deformed mini-Lpd + VASP condensates, `t = 5..180 s`, n >= 6). Calibrated from L-axis ticks (`x: 0..200 s`, `y: 0..1`), `res = (3.30, 1.25) px`. Note stray `"Time(s)"` label above panel B.
- **S3D**: Read — [`figs3d_frap_minilpd`, `figs3d_frap_actin`] (normalized FRAP intensity ± SD for mini-Lpd and actin in deformed mini-Lpd condensates, `t = 5..210 s`). Calibrated from L-axis ticks (`x: 0..250 s`, `y: 0..1`), `res = (1.40, 1.64) px`.
- **S3F**: Read — [`figs3f_frap_monomer_minilpd`, `figs3f_frap_actin`] (normalized FRAP intensity ± SD for monomer mini-Lpd and actin in deformed monomer mini-Lpd condensates, `t = 5..245 s`). Calibrated from L-axis ticks (`x: 0..250 s`, `y: 0..1`), `res = (1.49, 1.02) px`.
- **S3H**: Read — [`figs3h_frap_eps15`, `figs3h_frap_actin`] (normalized FRAP intensity ± SD for Eps15 [`t = 20..150 s`, pre-bleach `t = 5` occluded by red actin] and actin [`t = 5..150 s`] in deformed 1:1 WT Eps15:Eps15-Lifeact condensates). Calibrated from L-axis ticks (`x: 0..200 s`, `y: 0..1`), `res = (1.60, 1.23) px`.

### Figure S4 (p. 35) — *CytoSim model schematics*
- **S4A–J**: Left out — schematic drawings of CytoSim filaments, crosslinkers, and reaction rules (no quantitative data plots).

### Figure S5 (pp. 36–37) — *Viscosity and crosslinker copy-number sweeps*
- **S5A–D**: Left out — 3D CytoSim final snapshots (`t = 600 s` and `2400 s`) with no plotted axes.
- **S5E**: Read — [`figs5e_n20000_b0..b2`, `figs5e_n15000_b0..b2`, `figs5e_n10000_b0..b2`, `figs5e_n5000_b0..b2`, `figs5e_n1000_b0..b2`] (all 15 stacked-bar series of mini-Lpd dimer bound-filament fractions 0/1/2 ± SD across `N = 20000..1000` and `k_unbind = 0.001..10.0 1/s` at `k_bind = 1.0 1/s`, n = 5). Calibrated per box (`y: 0.0, 0.5, 1.0`), `ya_res <= 1.5 px`.

### Figure S6 (pp. 38–39) — *Mixed VASP:mini-Lpd & full 5×5 bivalent sweeps*
- **S6A, S6E**: Left out — 3D CytoSim final snapshots (`t = 600 s`).
- **S6B**: Read — [`figs6b_surf_sr_1000_0`, `figs6b_surf_sr_750_250`, `figs6b_surf_sr_500_500`, `figs6b_surf_sr_0_1000`, `figs6b_surf_rs_1000_0`, `figs6b_surf_rs_500_500`, `figs6b_surf_rs_250_750`, `figs6b_surf_rs_0_1000`, `figs6b_surf_rr_1000_0`, `figs6b_surf_rr_500_500`, `figs6b_surf_rr_0_1000`] (11 cleanly separable surface-area time series across Shell-Ring, Ring-Shell, Ring-Ring kinetics; `t = 550 s` omitted on `figs6b_surf_rs_1000_0` where red crosses blue). Left out curves: `sr 250:750` (occluded by `0:1000`), `rs 750:250` and `rr 750:250 / 250:750` (overlapping ring trajectories in `0.08–0.14`). Calibrated per box (`x: 0..600 s`, `y: 0.0..0.6`), `res <= (0.8, 1.4) px`.
- **S6C**: Read — [`figs6c_vasp_sr_b0..b4`, `figs6c_vasp_rs_b0..b4`, `figs6c_vasp_rr_b0..b4`] (all 15 stacked-bar series of VASP tetramer bound-filament fractions 0..4 ± SD across populated VASP counts `1000, 750, 500, 250`, n = 5). Calibrated per box on 6 major ticks (`0.0..1.0`), `ya_res <= 0.3 px`.
- **S6D**: Read — [`figs6d_lpd_sr_b0..b2`, `figs6d_lpd_rs_b0..b2`, `figs6d_lpd_rr_b0..b2`] (all 9 stacked-bar series of mini-Lpd dimer bound-filament fractions 0..2 ± SD across populated mini-Lpd counts `250, 500, 750, 1000`, n = 5; read from plotted green/orange/dark-red bars rather than the mis-colored legend, see `q_figs6d_miscolored_legend_and_missing_ylabels`). Calibrated per box on 6 major ticks (`0.0..1.0`), `ya_res <= 0.3 px`.
- **S6F**: Read — [`figs6f_surf_kb10_ku0p001`, `figs6f_surf_kb10_ku0p01`, `figs6f_surf_kb10_ku0p1`, `figs6f_surf_kb10_ku10`, `figs6f_surf_kb1_ku0p001`, `figs6f_surf_kb1_ku0p01`, `figs6f_surf_kb1_ku0p1`, `figs6f_surf_kb1_ku10`, `figs6f_surf_kb0p1_ku0p01`, `figs6f_surf_kb0p1_ku10`, `figs6f_surf_kb0p01_ku10`, `figs6f_surf_kb0p001_ku10`] (12 cleanly separable trajectories across the 5×5 `(k_bind, k_unbind)` grid; `t = 450 s` omitted on `figs6f_surf_kb10_ku0p01` and `t = 250 s` omitted on `figs6f_surf_kb1_ku0p01` where neighboring curves cross; overlapping curves — e.g. `ku = 1.0` light grey and all 5 identical shell trajectories in `k_bind = 0.01, 0.001` — are represented by the unobstructed top shell curve `ku = 10.0`). Calibrated per box (`x: 0..600 s`, `y: 0.0..0.6`), `res <= (1.1, 1.5) px`.
- **S6G**: Read — [`figs6g_dimers_kb10_b0..b2`, `figs6g_dimers_kb1_b0..b2`, `figs6g_dimers_kb0p1_b0..b2`, `figs6g_dimers_kb0p01_b0..b2`, `figs6g_dimers_kb0p001_b0..b2`] (all 15 stacked-bar series of bivalent mini-Lpd bound-filament fractions 0..2 ± SD across full 5×5 grid, n = 10; middle 3×3 agrees closely with Fig. 4K). Calibrated per box (`y: 0.0, 0.5, 1.0`), `ya_res <= 1.2 px`.

### Figure S7 (pp. 40–41) — *Dynamic multimerization model*
- **S7A**: Left out — 6-state monomer-configuration schematic.
- **S7B–C**: Left out — 3D CytoSim final snapshots (`t = 600 s`).
- **S7D**: Read — [`figs7d_surf_kf10_ks10`, `figs7d_surf_kf1_ks10`, `figs7d_surf_kf1_ks0p01`, `figs7d_surf_kf0p1_ks10`, `figs7d_surf_kf0p1_ks0p1`, `figs7d_surf_kf0p01_ks10`, `figs7d_surf_kf0p01_ks0p01`, `figs7d_surf_kf0p01_ks0p001`, `figs7d_surf_kf0p001_ks10`] (9 cleanly distinguishable trajectories across all 5 `k_multimer_form` rows, keyed by true physical `k_multimer_split` verified against Fig. S7C snapshots because the printed legend in Fig. S7D is inverted — see `q_figs7d_inverted_legend`; `t = 250 s` omitted on `figs7d_surf_kf1_ks10` before separation; overlapping curves inside the dense ring cluster at `kf = 10, 1.0, 0.1` and identical shell curves at `kf = 0.001` are left out where occluded). Calibrated per box (`x: 0..600 s`, `y: 0.0..0.6`), `res <= (1.0, 1.0) px`.
- **S7E**: Read — [`figs7e_kf10_fm0..bm2`, `figs7e_kf1_fm0..bm2`, `figs7e_kf0p1_fm0..bm2`, `figs7e_kf0p01_fm0..bm2`, `figs7e_kf0p001_fm0..bm2`] (all 30 stacked-bar series of the 6 monomer states `FM0, BM0, FM1, BM1, FM2, BM2` ± SD across the full 5×5 grid, n = 5). Calibrated per box from inner top/bottom spines (`1.0 / 0.0`, where 0.0 tick is omitted), `ya_res = 0.00 px`.
- **S7F**: Read — [`figs7f_multimer_len_kf10`, `figs7f_multimer_len_kf1`, `figs7f_multimer_len_kf0p1`, `figs7f_multimer_len_kf0p01`, `figs7f_multimer_len_kf0p001`] (median multimer chain length with actin filaments across all 5×5 `(k_multimer_form, k_multimer_split)` conditions from log10 violin long-dashed medians, with Q1/Q3 recorded in reading notes). Calibrated per box from long major log ticks (`10^2, 10^1, 10^0`), `ya_res <= 0.7 px`.
- **S7G**: Read — [`figs7g_multimer_len_kf10`, `figs7g_multimer_len_kf1`, `figs7g_multimer_len_kf0p1`, `figs7g_multimer_len_kf0p01`, `figs7g_multimer_len_kf0p001`] (median multimer chain length without actin filaments across all 5×5 conditions from log10 violin white median circles; all 25 medians resolve to exact integer chain lengths `1..5` within ±0.05). Calibrated per box from long major log ticks (`10^2, 10^1, 10^0`), `ya_res <= 0.9 px`.
