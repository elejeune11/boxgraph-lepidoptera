# Figure Inventory and Digitization Log — `freedman2019`

Complete page-by-page, panel-by-panel audit of every figure in `corpus/actin_bundling/2019_Freedman.pdf` (15 PDF pages: pp. 1–6 PNAS 116(33), 16192–16197 (2019); pp. 7–15 Supporting Information Appendix pp. SI-1–SI-9). All **39 quantitative figure readings** (**268 data points**) across **17 quantitative plot panels** were implemented as fast, self-contained `mode: code` scripts (`900–1200 dpi`) under `figures/` backed by shared calibration modules `figures/_calib.py` and `figures/_panels.py`. Every single reading was verified through BoxGraph (`boxgraph figure check` / `Builder.read_figure`), visually inspected across all 39 overlays stitched into compact contact sheets under `/tmp/figwork_freedman2019/` (all 268/268 markers, bars, error-bar caps, and curve points confirmed dead-center; 0 errors), and saved to `readings.json` (39/39 recipe hashes matching, 0 mismatches).

---

## Page 1 (printed p. 16192)
- No figures (Title, Abstract, Significance, and Introduction text only).

## Page 2 (printed p. 16193) — Figure 1: Cross-linker segregation in experiment and simulations
- **Fig. 1A**: Left out — experimental two-color TIRF fluorescence micrograph of F-actin bundled by fascin (red) and $\alpha$-actinin (cyan) (`2 um` scale bar); qualitative imaging with no quantitative axes.
- **Fig. 1B**: Left out — cartoon schematic of the coarse-grained AFINES bead-spring filament and cross-linker model; no plot axes.
- **Fig. 1C**: Left out — three qualitative AFINES bundle simulation snapshots at $\rho_{\text{short}}/\rho_{\text{long}} = 0.67, 1.0, 2.0$ (`2 um` scale bar); no quantitative axes.
- **Fig. 1D**: Left out — single-realization kymograph/barcode domain maps at $\rho_{\text{short}}/\rho_{\text{long}} = 0.50, 1.0, 2.0$ illustrating extracted short/long domains over `0..15 um` and `t = 0..2000 s`; ensemble-averaged domain lengths are digitized in Figs. 2, 3, 5, S2, and S5.

## Page 3 (printed p. 16194) — Figure 2: Effect of cross-linker and filament mechanics on domain length
- **Fig. 2A**: Left out — geometric schematic defining bending angle $\theta$ and contour gap length $l_g$ across a short/long cross-linker transition (Eq. 1); no quantitative axes.
- **Fig. 2B**: Left out — 1D equilibrium lattice Metropolis step schematic illustrating an initial `6 l_b` gap switching to `2 l_b` or `4 l_b` (Eq. 2; Boltzmann weight omission already logged in `q_fig2b_omits_initial_gap_energy`).
- **Fig. 2C (Mean domain length vs cross-linker length difference $\Delta l_{\text{xl}}$ at $L_p = 17\ \mu\text{m}$)**: Digitized both series with SEM error bars at `1200 dpi` (`cal_res = (0.22, 0.23) px`; log10 x-axis calibrated on ticks `[10, 30, 100, 300, 1000]` nm where $\Delta l_{\text{xl}} = 0\ \text{nm}$ is placed on the $\log_{10}(10\ \text{nm})$ tick custom-labeled `"0"`, linear y-axis `[0..6] um`):
  - `fig2c_lattice` (`fig2c_lattice_datum`, 7 black filled circles, `error_bars: sem`): `(0 nm, 0.825 ± 0.027)`, `(27, 1.16 ± 0.043)`, `(50, 1.39 ± 0.057)`, `(100, 2.02 ± 0.028)`, `(300, 2.98 ± 0.081)`, `(500, 3.43 ± 0.63)`, `(1000, 2.82 ± 0.87) um`.
  - `fig2c_afines` (`fig2c_afines_datum`, 7 red open circles, `error_bars: sem`): `(0 nm, 0.618 ± 0.11)`, `(27, 0.742 ± 0.12)`, `(50, 1.09 ± 0.14)`, `(100, 2.02 ± 0.49)`, `(300, 2.88 ± 0.72)`, `(500, 3.81 ± 1.1)`, `(1000, 4.51 ± 1.4) um`.
  - **Model comparisons (`equilibrium_lattice_mc_eq2`, `run_eq_lattice_fig2c_sweep`)**:
    - `cmp_fig2c_lattice_repro` (`reproduction`): **`not reproduced`** (3 `close`, 4 `not reproduced`). The Eq. 2 lattice MC increases monotonically (`0.69, 1.43, 1.61, 2.39, 4.30, 5.18, 6.19 um` at 200 sweeps; `6.64 um` at 400 sweeps), whereas plotted Fig. 2C peaks at `500 nm` (`3.43 um`) and drops non-monotonically to `2.82 ± 0.87 um` at `1000 nm` (logged in `q_fig2c_fig2d_lattice_discrepancy`).
    - `cmp_fig2c_lattice_vs_afines` (`approximation`): **`not reproduced`** (4 `reproduced`, 1 `close`, 2 `not reproduced` at `27` and `50 nm`).
- **Fig. 2D (Mean domain length vs filament persistence length $L_p$ at $\Delta l_{\text{xl}} = 100\ \text{nm}$)**: Digitized both series with SEM error bars at `1200 dpi` (`cal_res = (0.38, 0.20) px`; log10 x-axis ticks `[5, 10, 50, 100] um`, linear y-axis `[1.5..4.0] um`; all 7 markers sit exactly at $L_p = \{2.5, 7.5, 12.5, 17, 25, 75, 125\}\ \mu\text{m}$, matching $\kappa_B / 0.004$ and confirming `q_lp_kbt_rounding_movies3`):
  - `fig2d_lattice` (`fig2d_lattice_datum`, 7 black filled circles, `error_bars: sem`): `(2.5 um, 1.27 ± 0.021)`, `(7.5, 1.58 ± 0.030)`, `(12.5, 1.78 ± 0.020)`, `(17, 1.99 ± 0.042)`, `(25, 2.30 ± 0.052)`, `(75, 3.32 ± 0.14)`, `(125, 3.69 ± 0.32) um`.
  - `fig2d_afines` (`fig2d_afines_datum`, 7 red open circles, `error_bars: sem`): `(2.5 um, 2.11 ± 0.97)`, `(7.5, 1.94 ± 0.50)`, `(12.5, 1.89 ± 0.28)`, `(17, 2.41 ± 0.52)`, `(25, 2.56 ± 0.53)`, `(75, 2.98 ± 0.65)`, `(125, 3.02 ± 0.73) um`.
  - **Model comparisons (`equilibrium_lattice_mc_eq2`, `run_eq_lattice_fig2d_sweep`)**:
    - `cmp_fig2d_lattice_repro` (`reproduction`): **`not reproduced`** (2 `close`, 5 `not reproduced`) — Eq. 2 lattice MC predicts `1.55, 1.99, 2.35, 2.39, 2.81, 4.06, 4.50 um`, systematically ~20–32% above Fig. 2D's plotted black circles (`q_fig2c_fig2d_lattice_discrepancy`).
    - `cmp_fig2d_lattice_vs_afines` (`approximation`): **`close`** (4 `reproduced`, 3 `close`).

## Page 4 (printed p. 16195) — Figure 3 (KMC growing-bundle model) & Figure 4 (TIRF experiments)
- **Fig. 3A**: Left out — schematic of the 1D growing-bundle KMC model showing barbed-end site addition at rate $\mu_{\text{grow}}$ and cross-linker binding rates $\mu_{\text{on}}$ / $\mu_{\text{on}}^{\text{diff}}$; no plot axes.
- **Fig. 3B (Conditional dissimilar cross-linker switching probability vs gap distance $l_g$ in binding sites)**: Digitized all 20 filled blue circles (`n_g = 1..20`) at `1200 dpi` (`cal_res = (0.20, 0.00) px`):
  - `fig3b_switching_prob` (`fig3b_switching_prob_datum`, 20 points): `n_g = 1..3`: `-0.00015`, `n_g = 4`: `0.00454`, `n_g = 5`: `0.06333`, `n_g = 6`: `0.1746`, `n_g = 7`: `0.2735`, `n_g = 8`: `0.3421`, `n_g = 9`: `0.3874`, `n_g = 10`: `0.4174`, ..., `n_g = 20`: `0.4902`.
  - **Model comparison (`site_binding_probability_fig3b`, `run_fig3b_switching_prob_sweep`)**: `cmp_fig3b_switching_prob` (`reproduction`) -> **`not reproduced`** (5 `reproduced` at extreme tails `n=1..3, 19..20`, 4 `close` at `n=15..18`, **11 `not reproduced`** across `n=4..14`; e.g., at `n=6`, stated `(l_b = 37 nm, Delta l_xl = 27 nm)` predicts `0.2426` vs plotted `0.1746`).
  - **Forensic finding (`q_kmc_35nm_37nm_parameter_swap`)**: Across `n_g = 5..16`, the implied energy ratio $U_{\text{plot}}/U(37, 27)$ is constant at `1.363`, which equals $(37/35)^3 (29/27)^2 = 1.3629$ — proving Fig. 3B was computed with swapped `l_b = 35 nm` and `Delta l_xl = 37 - 8 = 29 nm`, which matches every unclipped point (`n_g = 5..16`) to `<= 0.0001`.
- **Fig. 3C (Left — Domain length vs $[\text{fascin}]/[\alpha\text{-actinin}]$ at $\mu_{\text{grow}} = 40\ \text{nm/s}$)**: Digitized all 3 series at `1200 dpi` (`cal_res = (0.17, 0.00) px`):
  - `fig3c_left_kmc_fascin` (`fig3c_left_kmc_fascin_datum`, 5 red filled circles): `(0.3, 0.787)`, `(0.5, 1.330)`, `(1.0, 3.999)`, `(1.5, 7.466)`, `(2.0, 10.03) um`.
  - `fig3c_left_kmc_alpha_actinin` (`fig3c_left_kmc_alpha_actinin_datum`, 5 cyan filled circles): `(0.3, 7.593)`, `(0.5, 6.039)`, `(1.0, 3.828)`, `(1.5, 3.154)`, `(2.0, 2.733) um`.
  - `fig3c_left_expt_fascin` (`fig3c_left_expt_fascin_datum`, 8 unoccluded/recoverable black squares from Ref. 16, `error_bars: unstated`): `(0.126, 0.54 ± 0.47)`, `(0.313, 2.39 ± 1.53)`, `(0.521, 1.10 ± 1.32)`, `(0.622, 2.63 ± 0.49)`, `(0.791, 1.85 ± 0.75)`, `(1.053, 5.80 ± 3.05)`, `(1.250, 3.48 ± 1.24)`, `(1.580, 6.22 ± 2.14) um`. (The 2 heavily occluded low-x squares near `0.18` and `0.26` are read cleanly on Fig. 5B's log10 axis; all 8 shared points agree with Fig. 5B to `<= 0.02 um`.)
  - **Model comparisons (`kmc_growing_bundle_conc_ratio`)**:
    - `cmp_fig3c_left_kmc_fascin` (`reproduction`): **`not reproduced`** (`0.78, 1.24, 3.24, 5.95, 8.57 um` vs plotted `0.79, 1.33, 4.00, 7.47, 10.03 um`; 15–20% low due to the `(35 nm, 29 nm)` swap).
    - `cmp_fig3c_left_kmc_alpha_actinin` (`reproduction`): **`not reproduced`** (`6.51, 4.93, 3.24, 2.52, 2.12 um` vs plotted `7.59, 6.04, 3.83, 3.15, 2.73 um`; 15–22% low).
    - `cmp_fig3c_left_kmc_vs_expt` (`validation`): **`close`** (6/8 `reproduced`, 2/8 `close`).
- **Fig. 3C (Right — Domain length vs polymerization rate $\mu_{\text{grow}}$ at equal $\mu_{\text{on}} = 0.2/\text{site/s}$)**: Digitized both series at `1200 dpi` (`cal_res = (0.00, 0.00) px`; note `mu_grow = 10 nm/s` is clipped above the top `12 um` axis frame and therefore omitted from `at`):
  - `fig3c_right_kmc_alpha_actinin` (`fig3c_right_kmc_alpha_actinin_datum`, 3 cyan circles): `(20 nm/s, 8.018)`, `(50, 3.140)`, `(100, 1.801) um`.
  - `fig3c_right_kmc_fascin` (`fig3c_right_kmc_fascin_datum`, 3 red circles peeking slightly behind cyan): `(20 nm/s, 8.145)`, `(50, 3.204)`, `(100, 1.828) um`.
  - **Model comparisons (`kmc_growing_bundle_fig3`, `run_kmc_fig3c_right_sweep`)**:
    - `cmp_fig3c_right_kmc_alpha_actinin` (`reproduction`): **`not reproduced`** (`6.32, 2.65, 1.64 um` vs plotted `8.02, 3.14, 1.80 um`; 9–21% low, again explained by `(35 nm, 29 nm)` giving `7.70, 3.24, 1.89 um`).
    - `cmp_fig3c_right_kmc_fascin` (`reproduction`): **`not reproduced`** (`6.39, 2.64, 1.66 um` vs plotted `8.15, 3.20, 1.83 um`).
- **Fig. 4A**: Left out — two-color TIRF fluorescence micrographs at `0.75 uM` (`40 min`) and `1.5 uM` (`20 min`) actin with `150 nM` fascin and `150 nM` $\alpha$-actinin (`2 um` scale bars); qualitative imaging.
- **Fig. 4B (Measured TIRF domain lengths at High vs Low labeled actin)**: Digitized both bar series with SEM error bars (`n = 2` replicates) at `900 dpi` (`cal_res = (0.00, 0.17) px`):
  - `fig4b_fascin` (`fig4b_fascin_datum`, 2 red bars, `error_bars: sem`): High Actin (`1.5 uM`, `46.99 nm/s`): `2.64 ± 0.20 um`; Low Actin (`0.75 uM`, `22.15 nm/s`): `1.35 ± 0.057 um`.
  - `fig4b_alpha_actinin` (`fig4b_alpha_actinin_datum`, 2 cyan bars, `error_bars: sem`): High Actin (`1.5 uM`): `2.21 ± 0.15 um`; Low Actin (`0.75 uM`): `8.29 ± 0.15 um`.

## Page 5 (printed p. 16196) — Figure 5: Competition between binding kinetics and polymerization in AFINES
- **Fig. 5A**: Left out — qualitative AFINES growing-bundle trajectory snapshots at varying $\rho_{\text{short}}$ (`2 um` scale bar); no plot axes.
- **Fig. 5B (Short cross-linker domain length vs density/concentration ratio at $k_{\text{grow}} = 40\ \text{nm/s}$)**: Digitized all 3 series at `1200 dpi` (`cal_res = (0.21, 0.22) px`; log10 x-axis `[0.1, 0.5, 1.0]`, linear y-axis `[0..10] um`):
  - `fig5b_afines_dl100` (`fig5b_afines_dl100_datum`, 7 red open circles, `error_bars: sem`): `(0.1, 0.270 ± 0.14)`, `(0.2, 0.553 ± 0.11)`, `(0.5, 0.895 ± 0.10)`, `(1.0, 1.89 ± 0.26)`, `(1.5, 4.36 ± 0.83)`, `(2.0, 5.62 ± 1.50)`, `(3.0, 8.31 ± 2.40) um`.
  - `fig5b_afines_dl27` (`fig5b_afines_dl27_datum`, 7 red filled circles, `error_bars: sem`): `(0.1, 0.274 ± 0.14)`, `(0.2, 0.308 ± 0.066)`, `(0.5, 0.520 ± 0.094)`, `(1.0, 0.904 ± 0.11)`, `(1.5, 1.52 ± 0.14)`, `(2.0, 2.37 ± 0.18)`, `(3.0, 5.00 ± 0.41) um`.
  - `fig5b_expt_fascin` (`fig5b_expt_fascin_datum`, all 10 black open squares from Ref. 16, `error_bars: unstated`): `(0.126, 0.563 ± 0.21)`, `(0.179, -0.008 ± 0)`, `(0.261, 0.590 ± 0.39)`, `(0.315, 2.40 ± 1.52)`, `(0.522, 1.10 ± 1.32)`, `(0.628, 2.61 ± 0.50)`, `(0.785, 1.86 ± 0.75)`, `(1.052, 5.81 ± 3.04)`, `(1.247, 3.48 ± 1.24)`, `(1.580, 6.22 ± 2.14) um`.
  - **Model validation (`kmc_growing_bundle_conc_ratio`, `run_kmc_fig5b_expt_sweep`)**: `cmp_fig5b_kmc_vs_expt` (`validation`) -> **`not reproduced`** (7/10 `reproduced`, 2/10 `close`, 1/10 `not reproduced` at the `x = 0.179` zero-length experimental point where KMC predicts `0.60 um`).
- **Fig. 5C (AFINES domain lengths vs $k_{\text{grow}}$ at $k_{\text{short}}^{\text{off},1}/k_{\text{long}}^{\text{off},1} = 1$)**: Digitized both series (9 points `k_grow = 2.5..100 nm/s`) with SEM error bars at `1200 dpi` (`cal_res = (0.05, 0.00) px`; vertical dashed experimental reference lines measure `22.15 nm/s` and `46.99 nm/s`, matching Fig. 4 caption to `0.01 nm/s`):
  - `fig5c_afines_short` (`fig5c_afines_short_datum`, 9 red open circles, `error_bars: sem`): `(2.5 nm/s, 2.38 ± 0.84)`, `(5, 2.83 ± 0.85)`, `(10, 1.82 ± 0.13)`, `(20, 1.84 ± 0.23)`, `(30, 1.89 ± 0.27)`, `(40, 2.22 ± 0.31)`, `(50, 1.92 ± 0.31)`, `(70, 1.88 ± 0.26)`, `(100, 1.88 ± 0.30) um`.
  - `fig5c_afines_long` (`fig5c_afines_long_datum`, 9 cyan open squares, `error_bars: sem`): `(2.5 nm/s, 3.88 ± 1.09)`, `(5, 2.84 ± 0.58)`, `(10, 1.96 ± 0.30)`, `(20, 1.41 ± 0.28)`, `(30, 1.16 ± 0.17)`, `(40, 1.06 ± 0.15)`, `(50, 1.02 ± 0.18)`, `(70, 0.917 ± 0.11)`, `(100, 0.886 ± 0.15) um`.
  - **Finding (`q_fig5c_crossover_rate_5nms_vs_3nms`)**: The curves cross directly at `k_grow = 5 nm/s` (`2.83` vs `2.84 um`), rather than text's `"~3 nm/s"`; both `2.5` and `5 nm/s` lie below Table S1's `l_g >= 0.008 um` (`8 nm/s`) bound.
- **Fig. 5D**: Left out — cartoon schematic defining cross-linker single-head ($k^{\text{on},1}, k^{\text{off},1}$) and double-head ($k^{\text{on},2}, k^{\text{off},2}$) transition rate constants (SI Eq. 4); no plot axes.
- **Fig. 5E (Top & Bottom — AFINES domain lengths vs $k_{\text{grow}}$ at single-filament unbinding ratios 2 and 4)**: Digitized all 4 series (`k_grow = 10..100 nm/s`) with SEM error bars at `1200 dpi` (`cal_res = (0.14, 0.00) px` on both stacked panels):
  - `fig5e_top_afines_short` (`fig5e_top_afines_short_datum`, 7 red open circles, ratio `2`): `(10 nm/s, 1.68 ± 0.30)`, `(20, 1.64 ± 0.23)`, `(30, 1.50 ± 0.19)`, `(40, 1.66 ± 0.24)`, `(50, 1.53 ± 0.16)`, `(70, 1.56 ± 0.16)`, `(100, 1.60 ± 0.21) um`.
  - `fig5e_top_afines_long` (`fig5e_top_afines_long_datum`, 7 cyan open squares, ratio `2`): `(10 nm/s, 3.15 ± 0.60)`, `(20, 1.92 ± 0.27)`, `(30, 1.59 ± 0.23)`, `(40, 1.38 ± 0.15)`, `(50, 1.32 ± 0.15)`, `(70, 1.19 ± 0.15)`, `(100, 1.08 ± 0.15) um`.
  - `fig5e_bot_afines_short` (`fig5e_bot_afines_short_datum`, 7 red open circles, ratio `4`): `(10 nm/s, 1.24 ± 0.26)`, `(20, 1.19 ± 0.18)`, `(30, 1.09 ± 0.18)`, `(40, 1.16 ± 0.17)`, `(50, 1.15 ± 0.15)`, `(70, 1.22 ± 0.18)`, `(100, 1.33 ± 0.22) um`.
  - `fig5e_bot_afines_long` (`fig5e_bot_afines_long_datum`, 7 cyan open squares, ratio `4`): `(10 nm/s, 4.38 ± 0.46)`, `(20, 2.83 ± 0.43)`, `(30, 1.96 ± 0.18)`, `(40, 1.94 ± 0.24)`, `(50, 1.78 ± 0.23)`, `(70, 1.43 ± 0.17)`, `(100, 1.13 ± 0.19) um`.

## Page 6 (printed p. 16197)
- No figures (Materials and Methods, Acknowledgments, and References only).

## Pages 7–10 (printed pp. SI-1–SI-4)
- No figures (SI cover page, SI Sections A–D equations/methods, and Table S1 parameter table on p. 10 / SI-4).

## Page 11 (printed p. SI-5) — SI Figures S1 & S2
- **Fig. S1**: Left out — illustrative single-trajectory, single-frame (`t = 1000 s`) 37-nm binning example showing how raw cross-linker coordinates along `0..15 um` are grouped into domains (SI Section B).
- **Fig. S2 (AFINES calibration — separate short and long domain lengths vs plotted length difference)**: Digitized both series with error bars (`error_bars: unstated`) at `900 dpi` (`cal_res = (0.53, 0.00) px`; log10 x-axis ticks `[0.01, 0.05, 0.10, 0.50, 1.0] um`, linear y-axis `[0..6] um`):
  - **Finding (`q_fig2c_figs2_log_axis_zero_offset`)**: Every point is plotted at $\Delta l_{\text{xl}} + 0.01\ \mu\text{m}$ (`0.01, 0.037, 0.06, 0.11, 0.31, 0.51, 1.01 um`) on a standard log10 axis without relabeling the ticks, shifting all 6 nonzero points rightward relative to their nominal values.
  - `figS2_afines_short` (`figS2_afines_short_datum`, 7 red open circles): `(0.01 um, 0.618 ± 0.072)`, `(0.037, 0.858 ± 0.067)`, `(0.06, 1.06 ± 0.076)`, `(0.11, 1.74 ± 0.22)`, `(0.31, 2.65 ± 0.40)`, `(0.51, 3.40 ± 0.52)`, `(1.01, 4.98 ± 0.88) um`.
  - `figS2_afines_long` (`figS2_afines_long_datum`, 7 cyan open squares): `(0.01 um, 0.617 ± 0.10)`, `(0.037, 0.631 ± 0.10)`, `(0.06, 1.14 ± 0.10)`, `(0.11, 2.31 ± 0.47)`, `(0.31, 3.12 ± 0.63)`, `(0.51, 4.20 ± 0.99)`, `(1.01, 4.03 ± 1.10) um`.

## Page 12 (printed p. SI-6) — SI Figure S3: Further details for the growing-bundle KMC model
- **Fig. S3A (Left — Waiting-time probability density $P(t)$ before same vs different cross-linker binding)**: Digitized both histogram curves at `1200 dpi` (`cal_res = (0.15, 0.15) px`):
  - `figS3a_left_same` (`figS3a_left_same_datum`, blue curve, 13 points `t = 0.5..8.0 s`): peaks at `(0.5 s, 0.2901 1/s)` and decays smoothly to `(8.0 s, 0.00524 1/s)`.
  - `figS3a_left_diff` (`figS3a_left_diff_datum`, orange curve, 15 points `t = 0.5..12.0 s`): rises from `(0.5 s, 0.1489 1/s)` to peak `(2.0 s, 0.2189 1/s)` and decays to `(12.0 s, 0.00382 1/s)`.
- **Fig. S3A (Right — Binding-site distance probability $P(l_g)$ for same vs different cross-linker binding)**: Digitized both histogram curves at `1200 dpi` (`cal_res = (0.20, 0.17) px`; printed inset values `N_diff/N_same = 0.016`, `epsilon = 4.1` already recorded in `data.yaml`):
  - `figS3a_right_same` (`figS3a_right_same_datum`, blue curve, 7 unoccluded sites `l_g = 4..10`; sites `1..3` rise above the `0.2` y-axis ceiling): `(4, 0.1340)`, `(5, 0.0743)`, `(6, 0.0375)`, `(7, 0.0167)`, `(8, 0.00747)`, `(9, 0.00435)`, `(10, 0.00300)`.
  - `figS3a_right_diff` (`figS3a_right_diff_datum`, orange curve, 14 sites `l_g = 2..18`): peaks around `l_g = 6–7` (`~0.0063`) and decays back toward zero by `l_g = 12..18`.
- **Fig. S3B (Bundle extension rate probability density $P(\mu_{\text{bund}})$)**: Digitized the cyan Type B curve (coincident with red Type A by symmetry) at `900 dpi` (`cal_res = (0.20, 0.00) px`; dashed vertical line sits at `40 nm/s`):
  - `figS3b_type_b` (`figS3b_type_b_datum`, 17 points `mu_bund = 5..245 nm/s`): rises from `(5 nm/s, 0.00522 s/nm)` to peak `(15–20 nm/s, 0.0170 s/nm)` and decays smoothly to `(245 nm/s, 0.00038 s/nm)`.
- **Fig. S3C (KMC master collapse $\epsilon_{\text{predicted}} = -\ln(N_{\text{diff}}/N_{\text{same}})$ vs $\mu_{\text{grow}}(\text{sites/s})/\mu_{\text{on}}(1/\text{s})$)**: Digitized all 7 $\mu_{\text{grow}}$ series (24 cleanly visible / crescent-recoverable markers; heavily occluded overlapping markers omitted rather than guessed) at `900 dpi` (`cal_res = (0.30, 0.20) px`; dashed Ref. 16 barrier line sits at `4.8`):
  - `figS3c_100nms` (`figS3c_100nms_datum`, 6 unoccluded pink circles): `(2.856, 5.435)`, `(5.710, 4.357)`, `(14.29, 2.615)`, `(28.57, 2.041)`, `(57.17, 1.423)`, `(142.8, 0.8895)`.
  - `figS3c_50nms` (`figS3c_50nms_datum`, 3 brown circles): `(1.429, 7.926)`, `(6.929, 3.640)`, `(71.33, 1.377)`.
  - `figS3c_40nms` (`figS3c_40nms_datum`, 4 purple circles): `(2.227, 5.805)`, `(5.719, 4.008)`, `(11.43, 2.826)`, `(22.84, 2.170)`.
  - `figS3c_30nms` (`figS3c_30nms_datum`, 6 red circles): `(0.858, 8.052)`, `(1.717, 7.200)`, `(4.290, 4.590)`, `(8.573, 3.292)`, `(16.84, 2.204)`, `(42.80, 1.623)`.
  - `figS3c_20nms` (`figS3c_20nms_datum`, 2 green circles): `(1.143, 7.301)`, `(2.862, 5.790)`.
  - `figS3c_10nms` (`figS3c_10nms_datum`, 1 orange crescent): `(1.427, 7.182)`.
  - `figS3c_5nms` (`figS3c_5nms_datum`, 2 blue circles): `(0.714, 8.119)`, `(1.427, 6.821)`.
  - **Model comparisons (`kmc_growing_bundle_rate_ratio`)**: All 7 comparisons (`cmp_figs3c_{100,50,40,30,20,10,5}nms`, `reproduction`) -> **`not reproduced`**: at the paper's stated `(l_b = 37 nm, Delta l_xl = 27 nm)`, computed $\epsilon_{\text{predicted}}$ sits systematically `0.1–0.6 k_B T` below the digitized master curve, while the measured x-coordinates (`2.856, 5.710, 14.29, 28.57, 57.17, 142.8`) equal `(100/35)/mu_on` to `< 0.06%` (`q_kmc_35nm_37nm_parameter_swap`).

## Page 13 (printed p. SI-7) — SI Figure S4: Additional TIRF data (unlabeled actin)
- **Fig. S4A**: Left out — wide-field two-color TIRF fluorescence micrographs of the same labeled-actin experiments as Fig. 4 (`4 um` scale bars); qualitative imaging.
- **Fig. S4B**: Left out — two-color TIRF fluorescence micrographs with labeled fascin and labeled $\alpha$-actinin but unlabeled actin (`1.5 uM` at `1200 s`, `3 uM` at `100 s`, `4 um` scale bars); qualitative imaging.
- **Fig. S4C (Measured TIRF domain lengths with unlabeled actin)**: Digitized both bar series with SEM error bars (`n = 2` replicates) at `900 dpi` (`cal_res = (0.00, 0.30) px`):
  - `figS4c_fascin` (`figS4c_fascin_datum`, 2 red bars, `error_bars: sem`): High Actin (`3 uM`): `1.83 ± 0.31 um`; Low Actin (`1.5 uM`): `1.81 ± 0.18 um`.
  - `figS4c_alpha_actinin` (`figS4c_alpha_actinin_datum`, 2 cyan bars, `error_bars: sem`): High Actin (`3 uM`): `1.62 ± 0.010 um`; Low Actin (`1.5 uM`): `2.48 ± 0.38 um`.

## Page 14 (printed p. SI-8) — SI Figure S5: AFINES fixed 1500–1600 s averaging window
- **Fig. S5 (AFINES short and long domain lengths vs $k_{\text{grow}}$ averaged over 1500–1600 s)**: Digitized both series (`k_grow = 0.01..0.10 um/s`) with SEM error bars at `900 dpi` (`cal_res = (0.02, 0.71) px`):
  - `figS5_afines_short` (`figS5_afines_short_datum`, 7 red open circles, `error_bars: sem`): `(0.01 um/s, 1.86 ± 0.15)`, `(0.02, 1.54 ± 0.21)`, `(0.03, 1.44 ± 0.19)`, `(0.04, 1.83 ± 0.26)`, `(0.05, 1.68 ± 0.20)`, `(0.07, 1.55 ± 0.21)`, `(0.10, 2.02 ± 0.39) um`.
  - `figS5_afines_long` (`figS5_afines_long_datum`, 7 cyan open squares, `error_bars: sem`): `(0.01 um/s, 1.91 ± 0.27)`, `(0.02, 1.83 ± 0.23)`, `(0.03, 1.54 ± 0.18)`, `(0.04, 1.79 ± 0.21)`, `(0.05, 1.67 ± 0.24)`, `(0.07, 2.05 ± 0.28)`, `(0.10, 2.26 ± 0.60) um`.

## Page 15 (printed p. SI-9)
- No figures (Movie S1–S7 captions and SI References 1–6 only).
