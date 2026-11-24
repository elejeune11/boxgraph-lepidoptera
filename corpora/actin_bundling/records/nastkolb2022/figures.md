# Figure digitization log — `nastkolb2022`

Source: T. Nast-Kolb, P. Bleicher, M. Payr & A. R. Bausch (2022), *"VASP localization to lipid bilayers induces polymerization driven actin bundle formation"*, Mol. Biol. Cell 33:ar91, 1–12 (`corpus/actin_bundling/2022_NastKolb.pdf`, PDF pages 1–18 including Supplemental Material pp. 13–18).

All 18 quantitative panels (9 main-text panels, 9 supplemental panels) were rendered at **1200 dpi** and digitized programmatically (`mode: code`, 58 readings in `figures.yaml`) with axis ticks calibrated directly from each panel image. Every reading was verified via `boxgraph figure check` and visually inspected on its overlay PNG in `figwork/<reading_id>.png`.

---

## Main-text figures (PDF pp. 3–8)

### Fig. 1 (PDF p. 3)
- **Fig. 1A**: *Left out* — schematic cartoon of VASP-mediated actin polymerization on an SLB (no quantitative data).
- **Fig. 1B**: *Left out* — TIRFM time-lapse micrographs (0–78 s at 50 nM VASP).
- **Fig. 1C**: *Left out* — TIRFM micrograph and single-bundle fluorescence kymograph (distance 0–10 µm vs time 100–300 s) with continuous intensity heatmap and no extracted quantitative trajectory.
- **Fig. 1D**: *Left out* — overlaid 0 s / 300 s TIRFM micrographs and temporal color-coded intensity-increase micrograph.

### Fig. 2 (PDF p. 4)
- **Fig. 2A, 2B, 2F**: *Left out* — TIRFM micrographs (0.125–0.5 µM VASP networks, 5–300 s swelling time series, and 1 µM VASP Atto488/Atto647N exchange).
- **Fig. 2C** (crop `[5900, 3520, 8880, 5580]`, 1200 dpi): *Read* (`fig2c_0s`, `fig2c_405s`).
  - **Calibration**: x-axis ticks 0–4 µm (5 unobstructed major ticks before legend border), linear misfit `rx = 1.20 px`; y-axis ticks 0–5×10⁶ AU (6 major ticks), linear misfit `ry = 2.34 px`.
  - **Curves read**: endpoint transverse profiles at `0 s` (darkest purple bottom trace; 20 positions 0.0–4.75 µm, using `1.2 µm` instead of `1.25 µm` where the 10 s purple trace overlaps 0 s) and `405 s` (brightest yellow top trace; 20 positions 0.0–4.75 µm), capturing the baseline vs swollen multi-peak bundle cross-section (e.g. central peak grows from ~5×10⁵–8×10⁵ AU at 0 s to ~4.3×10⁶ AU across 2.0–3.0 µm at 405 s).
  - **Left out within panel**: the 11 intermediate time traces (10–300 s) share a continuous viridis color ramp and overlap densely on the steep bundle flanks, preventing unambiguous per-trace separation.
- **Fig. 2D** (crop `[944, 5504, 4984, 8344]`, 1200 dpi): *Read* (8 readings: histogram markers `fig2d_hist_005`, `fig2d_hist_0125`, `fig2d_hist_025`, `fig2d_hist_050`; log-normal fit curves `fig2d_fit_005`, `fig2d_fit_0125`, `fig2d_fit_025`, `fig2d_fit_050`).
  - **Calibration**: x-axis ticks 0–5 µm (6 major ticks), `rx = 1.16 px`; y-axis ticks 0.0–0.6 (7 major ticks), `ry = 1.29 px`.
  - **Notes**: 0.2 µm bin centers. For `0.05 µM VASP` (`fig2d_hist_005`), bins `0.2–0.8 µm` are read cleanly (`0.2 µm` partially occluded by `0.125 µM` red circle, interval widened accordingly); bins `1.0` and `1.2 µm` are occluded behind the `0.125 µM` red circles and omitted. For `fig2d_fit_0125`, `x = 0.5 µm` is covered by the `0.05 µM` dark gray curve, so unobstructed positions `0.25, 0.3, 0.35, 0.45, 0.65, 0.7, 0.9, 1.1 µm` are read.
- **Fig. 2E** (crop `[4948, 5504, 8828, 8344]`, 1200 dpi): *Read* (`fig2e_bundle_width`, `fig2e_bundle_width_median`, `fig2e_pore_size`).
  - **Calibration**: x-axis ticks 0.00–1.00 µM (5 major ticks); left y-axis 0–8 µm (9 major ticks, `ry = 1.07 px`); broken right y-axis lower segment 0–7 µm² (8 major ticks) + upper segment 17–18 µm² with 17.5 minor tick (`ry = 1.45 px`, bridging axis break cleanly).
  - **Key readings**:
    - Bundle width medians [IQR; 1–99th whiskers] (µm): 0.05 µM: `0.547` (`[0.513, 0.581]`; `[0.259, 1.034]`, thin single-band box so median uncertainty spans full IQR); 0.125 µM: `0.581` (`[0.513, 0.636]`; `[0.259, 1.034]`); 0.25 µM: `0.945 ± 0.012` (`[0.581, 1.328]`; `[0.259, 2.473]`, closely reproducing printed `0.942 µm`); 0.50 µM: `1.589 ± 0.012` (`[1.034, 2.350]`; `[0.259, 4.234]`); 1.00 µM: `4.357 ± 0.012` (`[3.254, 5.358]`; `[0.787, 7.760]`).
    - Pore size medians [1–99th whiskers] (µm²): 0.05 µM: `1.393` (`[0.015, 17.16]`); 0.125 µM: `0.998` (`[0.015, 7.115]`); 0.25 µM: `0.594` (`[0.015, 3.929]`); 0.50 µM: `0.436` (`[0.015, 2.147]`); 1.00 µM: `0.199` (`[0.015, 1.963]`).
  - **Dual bundle-width readings**: `fig2e_bundle_width` attaches the 1–99th percentile population whiskers (`error_bars: range`), whereas `fig2e_bundle_width_median` preserves tight reading precision around the medians to test the paper's linear scaling claim `W ~ c_VASP` (`bundle_width_vasp_scaling`).

### Fig. 3 (PDF p. 5)
- **Fig. 3A–3D**: *Left out* — schematics and TIRFM micrographs of VASP_ΔTetra, α-actinin1b, mDia1, and mDia1 + α-actinin1b.
- **Fig. 3E** (crop `[6440, 2920, 9120, 5000]`, 1200 dpi): *Read* (all 6 series: `fig3e_vasp`, `fig3e_vasp_dtetra`, `fig3e_actinin_050`, `fig3e_actinin_005_nocp`, `fig3e_mdia1`, `fig3e_mdia1_actinin`).
  - **Calibration**: x-axis ticks 0–20 min (11 major ticks), `rx = 1.50 px`; y-axis ticks 0.0–1.0 (6 major ticks), `ry = 0.40 px`.
  - **Key finding**: `fig3e_mdia1` (all 21 unobstructed brown hexagons, 1.08–18.87 min) peaks at `0.95` around `6.88–8.88 min` and declines thereafter to `0.66` at `18.87 min`, confirming `q_figure_panel_misreferences` (p. 6 cites "after 6 min (Figure 3D)" instead of Fig. 3E, where decline begins after ~8–9 min).
- **Fig. 3F** (crop `[6440, 4944, 9120, 7084]`, 1200 dpi): *Read* (`fig3f_bundle_width`, `fig3f_pore_size`).
  - **Calibration**: 5 categorical conditions (`1 = VASP`, `2 = VASP_ΔTetra`, `3 = α-actinin`, `4 = mDia1`, `5 = mDia1 + α-actinin`); left y-axis 0–6 µm (7 major ticks, `ry = 1.06 px`); broken right y-axis lower 0–20 µm² (11 major ticks) + upper 44–46 µm² using tick 46 and minor tick 45 (`ry = 1.56 px`).
  - **Key readings & settled question**:
    - Bundle width medians [1–99th whiskers] (µm): VASP `0.938` (`[0.262, 2.476]`, matches printed `0.942 µm`), VASP_ΔTetra `0.588` (`[0.186, 3.092]`, matches printed `0.584 µm`), α-actinin `0.523` (`[0.264, 0.833]`), mDia1 `0.773` (`[0.252, 1.453]`), mDia1+α-actinin `0.773` (`[0.252, 4.925]`). Directly confirms `q_max_bundle_width_3um_vs_fig3f_whisker`: tetrameric VASP's 99th-percentile whisker in Fig. 3F is `2.48 µm`, not `3 µm` (only the Fig. 2D histogram tail reaches `3.0 µm`), whereas VASP_ΔTetra reaches `3.09 µm` (~3.2 µm printed).
    - Pore size medians [1–99th whiskers] (µm²): VASP `0.580` (`[0.099, 3.901]`), VASP_ΔTetra `2.799` (`[0.00, 21.13]`), α-actinin `0.454` (`[0.112, 3.594]`), mDia1 `0.462` (`[0.013, 45.06]`), mDia1+α-actinin `0.210` (`[0.00, 13.12]`).

### Fig. 4 (PDF p. 7)
- **Fig. 4A, 4D**: *Left out* — TIRFM micrographs of VASP/actin colocalization on SLB/ML and pre-actin VASP phase separation.
- **Fig. 4B** (crop `[864, 4504, 4904, 7484]`, 1200 dpi): *Read* (all 4 series: `fig4b_vasp_slb`, `fig4b_vasp_ml`, `fig4b_dtetra_slb`, `fig4b_dtetra_ml`).
  - **Calibration**: x-axis ticks 0–40 min (5 major ticks), `rx = 0.00 px`; y-axis ticks 0.0–1.0 (6 major ticks), `ry = 2.41 px`.
  - **Notes**: VASP on SLB (11 squares, 0.03–37.98 min) rises from `0.623` to plateau `0.95–0.96`; VASP_ΔTetra on SLB (6 diamonds, 0–41.99 min) rises from `0.411` to `0.83–0.88`; both ML series remain near `0.05–0.12` (note `t = 0` red triangle of VASP on ML is occluded behind the green VASP_ΔTetra ML triangle, so 6 unobstructed points 1.99–44.97 min are read).
- **Fig. 4C** (crop `[4920, 4520, 8960, 7500]`, 1200 dpi): *Read* (all 4 curves: `fig4c_slb_prepol`, `fig4c_slb_pol30`, `fig4c_ml_prepol`, `fig4c_ml_pol30`).
  - **Calibration**: x-axis fitted on 7 unobstructed major ticks 0.0–1.2 min (`rx = 0.00 px`); y-axis 0.0–1.2 (7 major ticks, `ry = 0.00 px`); sampled at 15 post-bleach times `0.05–1.30 min` (pre-bleach flat baseline `-0.2–0.0 min` at 1.0 and vertical bleach drop at `t = 0` excluded so sweeps test the post-bleach exponential recovery model).
  - **Key endpoints at `t = 1.30 min`**: SLB prePol `0.962 ± 0.012` (printed `0.97`), SLB 30 min Pol `0.641 ± 0.006` (printed `0.65`), ML prePol `0.607 ± 0.015` (printed `0.61`), ML 30 min Pol `0.535 ± 0.007` (printed `0.57`). Directly resolves `q_frap_amplitude_vs_mobile_fraction_reading`: on ML prePol, `0.61` matches the `t = 1.3 min` window endpoint (`0.607`), whereas `0.61*(1 - 2^(-t/0.399))` reaches only `0.546` at `1.3 min` and systematically underpredicts the entire curve (an asymptote `a ≈ 0.68` is needed with `t_1/2 = 0.399 min`).

### Fig. 5 (PDF p. 8)
- **Fig. 5A, 5D**: *Left out* — TIRFM micrographs of actin/Texas Red-DHPE lipids and mechanistic schematic.
- **Fig. 5B** (crop `[4000, 872, 6410, 2912]`, 1200 dpi): *Read* (`fig5b_lipids`, `fig5b_actin`).
  - **Calibration**: x-axis 0–8 µm (5 major ticks), `rx = 1.13 px`; left y-axis 0.90–1.30 (8 ticks present in PDF—tick `1.15` is missing a tick line in the source PDF, handled explicitly—`ry = 1.12 px`); right y-axis 0.60–1.00 (9 major ticks, `ry = 1.13 px`). Sampled at 34 positions 0.00–8.25 µm showing inverse correlation between actin peaks (up to `0.985` at 7.0 µm) and lipid troughs (down to `0.908` at 7.0 µm).
- **Fig. 5C** (crop `[6420, 880, 9120, 3000]`, 1200 dpi): *Read* (`fig5c_actin`, `fig5c_vasp`, `fig5c_lipids`), all with `error_bars: sd` (n = 5 bundles, p. 11).
  - **Calibration**: x-axis 0–60 min (7 major ticks), `rx = 0.21 px`; left red axis `-0.4` to `1.0` (`ry = 0.95 px`), inner right green axis `1.3` to `2.2` (`ry = 1.07 px`), outer right blue axis `0.90` to `1.25` (`ry = 0.95 px`).
  - **Key readings & settled question**:
    - Actin (`fig5c_actin`, 14 SD bars 1.54–60.60 min; `t = 0` point sits on the left spine without an SD bar): rises from `0.192 ± 0.031` (1.54 min) to plateau `0.945 ± 0.025` (25.3 min) and `0.867 ± 0.042` (60.6 min).
    - VASP ratio (`fig5c_vasp`, 14 SD bars 0.03–60.79 min): rises from `1.47 ± 0.25` (0.03 min) to peak `1.95 ± 0.17` (5.11 min), relaxing to `1.62 ± 0.11` (60.8 min).
    - Lipid ratio (`fig5c_lipids`, 14 SD bars 0–60.6 min): drops from `0.977 ± 0.029` (t ≈ 0) to minimum `0.910 ± 0.025` (4.81 min, ~9% depletion below 1.00) and relaxes to `0.951 ± 0.017` (60.6 min, ~5% below 1.00). Quantitatively confirms `q_lipid_relaxation_half_intensity_vs_fig5c`: lipid fluorescence never drops to 50% of initial intensity; rather, ~half of the ~9% initial depletion relaxes by 60 min.

---

## Supplemental figures (PDF pp. 15–18; captions on p. 14)

### Fig. S1 (PDF p. 15)
- **Fig. S1A–S1F**: *Left out* — control TIRFM micrographs (no VASP SLB/ML, no Ni-NTA, bulk 0.5 µM VASP, Ni-agarose beads) and confocal XZ cross-section micrograph.

### Fig. S2 (PDF p. 16)
- **Fig. S2A–S2E**: *Left out* — TIRFM micrographs across 0.05–1 µM VASP.
- **Fig. S2F** (crop `[5340, 2720, 8014, 5100]`, 1200 dpi): *Read* (all 5 curves: `figs2f_005`, `figs2f_0125`, `figs2f_025`, `figs2f_050`, `figs2f_100`; 14 times 0.5–60 min).
  - **Calibration**: x-axis 0–60 min (7 major ticks), `rx = 0.21 px`; y-axis 0–10000 AU (6 major ticks), `ry = 0.19 px`. Plateau intensities at 60 min: `562, 738, 6557, 7803, 9432 AU` for 0.05, 0.125, 0.25, 0.5, 1 µM VASP.
- **Fig. S2G** (crop `[0, 5264, 2660, 7440]`, 1200 dpi): *Read* (`figs2g_decoration`; `error_bars: sd`, p. 11).
  - **Calibration**: x-axis 0.0–1.0 µM (6 major ticks), `rx = 0.20 px`; y-axis 0.0–1.0 (11 major ticks), `ry = 0.32 px`.
  - **Values (mean ± SD)**: 0.05 µM: `0.376 ± 0.037`; 0.125 µM: `0.454` (error bar hidden inside marker, bounded by marker extent `[0.441, 0.466]`); 0.25 µM: `0.651 ± 0.049`; 0.50 µM: `0.760 ± 0.072`; 1.00 µM: `0.950 ± 0.038`.
- **Fig. S2H** (crop `[2660, 5280, 5340, 7456]`, 1200 dpi): *Read* (`figs2h_vasp_fluor`; `error_bars: unstated`).
  - **Calibration**: x-axis 0.0–2.0 µM (5 major ticks), `rx = 0.30 px`; y-axis 0–600 AU (7 major ticks), `ry = 0.21 px`.
  - **Values**: 0.05 µM: `4.7` (bar `[1.0, 9.6]`); 0.125 µM: `47.4` (`[24.6, 70.3]`); 0.25 µM: `200` (`[172, 228]`); 0.50 µM: `278` (`[257, 299]`); 1.00 µM: `380` (`[334, 426]`); 2.00 µM: `507` (`[490, 524]`).
- **Fig. S2I** (crop `[5288, 5188, 8002, 7444]`, 1200 dpi): *Read* (`figs2i_hist_100`, `figs2i_fit_100`).
  - **Calibration**: x-axis 0–8 µm (5 major ticks), `rx = 0.15 px`; y-axis 0.0–0.6 (7 major ticks), `ry = 0.43 px`.
  - **Notes**: Only the new `1 µM VASP` purple diamond histogram (42 bins, 0.2–8.4 µm) and purple log-normal fit curve (16 points, 0.5–8.1 µm; read after exact L1 diamond subtraction to prevent diamond-vertex leakage) are digitized here, since the `0.05–0.5 µM` series replot the exact Fig. 2D data at lower x-resolution.

### Fig. S3 (PDF p. 17)
- **Fig. S3A–S3D**: *Left out* — TIRFM micrographs at 0.5, 1, 2, 4 µM actin (0.25 µM VASP).
- **Fig. S3E** (crop `[2660, 2680, 5360, 5182]`, 1200 dpi): *Read* (all 4 series: `figs3e_05`, `figs3e_1`, `figs3e_2`, `figs3e_4`).
  - **Calibration**: x-axis 0–60 min (7 major ticks), `rx = 0.00 px`; y-axis 0–6000 AU (7 major ticks), `ry = 0.00 px`. Final intensities near 58.5 min: `1502, 2904, 4793, 5249 AU` for 0.5, 1, 2, 4 µM actin (0.5 and 1 µM series start at `2 min` because `1 min` markers are occluded by 2 µM/4 µM traces and slightly negative from background subtraction).
- **Fig. S3F** (crop `[5360, 2664, 8037, 5166]`, 1200 dpi): *Read* (`figs3f_bundle_width`, `figs3f_pore_size`).
  - **Calibration**: left y-axis 0.0–3.0 µm (7 major ticks, `ry = 0.23 px`); right y-axis 0–12 µm² (7 major ticks, `ry = 0.32 px`).
  - **Values**: Bundle width medians [1–99th] (µm): 0.5 µM `0.574` (`[0.256, 1.073]`), 1 µM `0.770` (`[0.256, 1.561]`), 2 µM `0.942` (`[0.261, 2.479]`, matching printed `0.942 µm` to 0.0001 µm!), 4 µM `0.812` (`[0.256, 1.643]`). Pore size medians [1–99th] (µm²): 0.5 µM `2.109` (`[0.033, 11.18]`), 1 µM `1.229` (`[0.016, 6.224]`), 2 µM `0.599` (`[0.016, 3.924]`), 4 µM `0.382` (`[0.016, 1.502]`).

### Fig. S4 (PDF p. 18)
- **Fig. S4A, S4B, S4D, S4E, S4G, S4I**: *Left out* — TIRFM micrographs (SLB/ML without CP, prepolymerized F-actin on SLB/ML, and pre-actin VASP clusters).
- **Fig. S4C** (crop `[5309, 108, 8038, 2648]`, 1200 dpi): *Read* (6 readings: histograms `figs4c_hist_vasp`, `figs4c_hist_nocp_slb`, `figs4c_hist_nocp_ml`; log-normal fits `figs4c_fit_vasp`, `figs4c_fit_nocp_slb`, `figs4c_fit_nocp_ml`).
  - **Calibration**: x-axis 0.0–3.0 µm (7 major ticks), `rx = 0.21 px`; y-axis 0.0–0.5 (6 major ticks), `ry = 0.40 px`.
- **Fig. S4F** (crop `[5320, 2720, 8049, 5180]`, 1200 dpi): *Read* (`figs4f_bundle_width`, `figs4f_pore_size`).
  - **Calibration**: 3 categorical conditions (`1 = VASP`, `2 = VASP noCP SLB`, `3 = VASP noCP ML`); left y-axis 0.0–3.0 µm (`ry = 0.14 px`); broken right y-axis lower 0–4 µm² + upper 18.5–19 µm² (`ry = 0.30 px`).
  - **Values**: Bundle width medians [1–99th] (µm): VASP `0.942` (`[0.261, 2.478]`), VASP noCP SLB `0.931` (`[0.258, 1.966]`), VASP noCP ML `0.577` (`[0.258, 1.291]`). Pore size medians [1–99th] (µm²): VASP `0.597` (`[0.016, 3.923]`), VASP noCP SLB `0.416` (`[0.016, 2.016]`), VASP noCP ML `0.500` (`[0.016, 18.24]`).
- **Fig. S4H** (crop `[0, 7560, 3000, 10044]`, 1200 dpi): *Read* (`figs4h_cluster_area`; `error_bars: unstated`).
  - **Calibration**: x-axis 0.0–2.0 µM (5 major ticks), `rx = 0.20 px`; y-axis 0.0–1.0 (6 major ticks), `ry = 0.30 px`.
  - **Values**: 0.05 µM: `0.007` (bar hidden within marker, `[0.00, 0.022]`); 0.125 µM: `0.044` (bar hidden within marker, `[0.029, 0.058]`); 0.25 µM: `0.176` (`[0.103, 0.250]`); 0.50 µM: `0.431` (`[0.258, 0.604]`); 1.00 µM: `0.596` (`[0.432, 0.761]`); 2.00 µM: `0.591` (`[0.308, 0.874]`). Quantifies `q_vasp_phase_separation_threshold_tension`: nonzero cluster area fraction already reaches `0.176` at `0.25 µM` despite main text p. 6 stating `>= 0.5 µM`.

---

## Model sweeps & comparison verdicts against digitized figures

| Comparison ID | Purpose | Model vs. Digitized Datum | Verdict | Key quantitative finding |
| --- | --- | --- | --- | --- |
| `cmp_width_scaling_fig2e_median` | `validation` | `bundle_width_vasp_scaling` vs. `bundle_width_median_vs_vasp_fig2e` | **not reproduced** (0/5 reproduced, 1 close) | Best-fit proportionality constant `4.843 µm/µM` predicts `[0.242, 0.605, 1.211, 2.421, 4.843] µm` vs. measured medians `[0.547, 0.581, 0.945, 1.589, 4.357] µm`; fails due to the ~0.55–0.58 µm diffraction/segmentation floor at 0.05–0.125 µM and superlinear jump at 1.0 µM (95% bilayer coverage). |
| `cmp_width_scaling_fig2e_range` | `validation` | `bundle_width_vasp_scaling` vs. `bundle_width_vs_vasp_fig2e` | **close** (4/5 reproduced, 1 close) | Inside the wide 1–99th percentile whiskers at 0.125–1.00 µM, but `0.242 µm` at 0.05 µM still lies slightly below the 1st-percentile lower whisker (`0.259 µm`). |
| `cmp_frap_fit_slb_prepol_fig4c` | `fit_quality` | `frap_recovery_fit` vs. `frap_curve_slb_prepol_fig4c` | **not reproduced** (3/15 reproduced, 4 close) | Using printed `a = 0.97`, `t_1/2 = 0.102 min`, `c = 0`, residuals show classic 2D-diffusion vs. single-exponential curvature: underpredicts early rise (`-0.048` at 0.05 min), overshoots mid-recovery (`+0.069` at 0.4 min), and converges at late plateau (`0.970` vs. `0.962` at 1.3 min). |
| `cmp_frap_by_thalf_slb_prepol_fig4c` | `fit_quality` | `frap_recovery_by_half_time` vs. `frap_curve_slb_prepol_fig4c` | **not reproduced** (3/15 reproduced, 4 close) | Identical curve and verdict to `cmp_frap_fit_slb_prepol_fig4c`, propagating printed rounding of both `0.97` and `0.102 min`. |
| `cmp_frap_fit_ml_prepol_fig4c` | `fit_quality` | `frap_recovery_fit` vs. `frap_curve_ml_prepol_fig4c` | **not reproduced** (1/15 reproduced, 3 close) | Using printed `a = 0.61`, `t_1/2 = 0.399 min`, `c = 0`, computed curve underpredicts the measured ML prePol trace at all 15 time points (`-0.006` to `-0.073`), reaching `0.546` vs. `0.607 ± 0.015` at 1.3 min—confirming printed `0.61` is the 1.3 min window endpoint rather than the exponential asymptote `a` (~0.68). |
| `cmp_frap_by_thalf_ml_prepol_fig4c` | `fit_quality` | `frap_recovery_by_half_time` vs. `frap_curve_ml_prepol_fig4c` | **not reproduced** (1/15 reproduced, 3 close) | Identical curve and verdict to `cmp_frap_fit_ml_prepol_fig4c`. |

