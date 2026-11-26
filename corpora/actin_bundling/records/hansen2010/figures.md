# Figure Digitization Log — Hansen & Mullins (2010)

Source PDF: `2010_Hansen.pdf` (25 pages: main article pp. 1–14 / printed pp. 571–584; Supplemental Material Figs. S1–S5 & Video legends pp. 15–25 / printed pp. S1–S5 and Fig. S5 pp. 1–5).

- **Summary:** Digitized **84 figure readings (698 data points)** spanning all **32 quantitative plot panels** with calibrated numeric axes across main text Figs. 1–6 and Supplemental Figs. S1–S5, all rendered at **600 dpi** (2× native 300 ppi main/Fig. S5 JPEGs; 1× native 600 ppi Figs. S1–S4 JPEGs). Every reading runs self-contained image-based frame/tick calibration and detection in `figures/*.py` + `figures/_calib.py`, passed `boxgraph figure check` with 0 errors, and was visually verified point-by-point on compact overlay contact sheets in `figwork/`. Saved frozen readings in `readings.json`.

---

## Main Text Figures (Figs. 1–7)

### Figure 1 (PDF p. 3 / printed p. 573)
- **Fig. 1A:** Left out — domain-architecture cartoon schematic of Cy3-KCK-hVASP (no numeric axes).
- **Fig. 1B:** Left out — SYPRO Red–stained SDS-PAGE cosedimentation gel (qualitative gel bands; subunit MWs already recorded as scalar data).
- **Fig. 1C (READ — 3 readings, 48 pts):** `fig1c_actin_alone` (16 pts), `fig1c_vasp` (16 pts), `fig1c_kck_vasp` (16 pts). Pyrene-actin polymerization time courses (0–1000 s, fluorescence 0–3500 a.u.) calibrated from inward frame ticks (`< 0.5 px` residual). Traced blue, black, and red curves cleanly after masking the `"Actin alone"` pointer line; confirms identical acceleration by unlabeled VASP and KCK-VASP (63 nM tetramer).
- **Fig. 1D (READ — 3 readings, 14 pts + SD bars):** `fig1d_50mM_vasp` (6 pts, 0–100 nM), `fig1d_50mM_cy3_vasp` (5 pts, 1–100 nM), `fig1d_100mM_vasp` (3 pts, 0, 50, 100 nM; no Cy3-VASP bars at 100 mM KCl). Bar tops and SD error-bar caps verified visually. Both 50 mM KCl series validate the Fig. 4G Michaelis-Menten model (`cmp_fig1d_50mM_vasp_rates` = `reproduced`, `cmp_fig1d_50mM_cy3_vasp_rates` = `reproduced`). At 100 mM KCl, rates rise `8.29 -> 9.42 -> 10.28 sub/s` (`+13.5%` at 50 nM, `+23.9%` at 100 nM), confirming `q_audit_acceleration_100mM_quantity`.

### Figure 2 (PDF p. 4 / printed p. 574)
- **Fig. 2A:** Left out — TIRF fluorescence micrographs of Alexa 488 F-actin and Cy3-VASP max-intensity projection.
- **Fig. 2B:** Left out — single representative 4-step photobleaching fluorescence trace (arbitrary fluorescence units; single-Cy3 bleaching `tau_1 = 11.5 ± 0.07 s` is already recorded as a scalar datum).
- **Fig. 2C:** Left out — TIRF kymograph of 0.25 nM Cy3-VASP interacting with F-actin.
- **Fig. 2D (READ — 2 readings, 21 pts):** `fig2d_survival_fit` (12 pts, 0.5–50 s, red double-exponential fit curve) and `fig2d_survival_data` (9 pts, 0.5–15 s, black empirical 1-CF data curve). Calibrated from inward spine ticks. `cmp_fig2d_survival_data_curve` = `close`; on `cmp_fig2d_survival_fit_curve`, `t = 3–50 s` matches the printed `86%/14%, tau_1 = 0.36 s, tau_2 = 4.99 s` curve to `0.0007–0.013`, while `t = 1 s` sits slightly above due to steep near-spine slope (`-1.6 / s` over ~20 px).
- **Fig. 2E & inset (READ — 2 readings, 42 pts):** `fig2e_diffusive_hist` (30 histogram bars, 0.05–1.50 s, counts 0–200) and `fig2e_inset_1cf_fit` (12 pts, 0.1–5.0 s, red single-exponential 1-CF fit curve). Compared against `tau_1 = 0.44 s` in `cmp_fig2e_inset_1cf_curve` (10/12 points agree within pixel interval).
- **Fig. 2F:** Left out — TIRF fluorescence micrographs comparing tetramer, monomer, dimer, and GAB/FAB mutants on F-actin.
- **Fig. 2G:** Left out — TIRF kymographs of wild-type and mutant Cy3-VASP on F-actin.
- **Fig. 2H (READ — 4 readings, 34 pts):** `fig2h_vasp_latb` (10 pts, black), `fig2h_vasp_lil_aaa` (7 pts, red dashed), `fig2h_vasp_rk_ee` (8 pts, blue dotted), `fig2h_vasp_2uM_actin` (9 pts, green dashed). **Key finding (`q_fig2h_vasp_latb_curve_mismatch`):** `VASP^(RK->EE)` (`cmp_fig2h_vasp_rk_ee_curve` = `reproduced`) and `+2 uM Actin` (`cmp_fig2h_vasp_2uM_actin_curve` = `reproduced`) reproduce their Fig. 2H inset table rows (and `LIL->AAA` agrees for `t >= 0.3 s`), whereas the black `VASP +LatB` curve sits `0.032–0.071` below its Fig. 2H table row (`82%, 0.41 s, 5.19 s`; `cmp_fig2h_vasp_latb_vs_table` = `not reproduced`) and instead reproduces the Fig. 2D no-LatB curve (`86%, 0.36 s, 4.99 s`; `cmp_fig2h_vasp_latb_vs_fig2d` = `reproduced` at all 10 points).

### Figure 3 (PDF p. 5 / printed p. 575)
- **Fig. 3A:** Left out — TIRF kymographs of barbed-end capture ± 2 uM actin/LatB.
- **Fig. 3B:** Left out — dual-color Alexa 488/568 F-actin + Cy3-VASP TIRF micrographs.
- **Fig. 3C (READ — 4 readings, 40 pts):** `fig3c_v`, `fig3c_v_cp`, `fig3c_v_a`, `fig3c_v_a_cp` (10 filament pixel positions each, 1–10; probability 0–0.25). Captures the sharp pixel-1 barbed-end peak (~0.244) in `+Actin` vs uniform ~0.09–0.11 localization in VASP alone or +5 nM CP.
- **Fig. 3D (READ — 2 readings, 4 pts + SEM bars):** `fig3d_barbed_end` (2 bars: `+Actin`, `+Actin, CP`) and `fig3d_lateral_side` (2 bars: `+Actin`, `+Actin, CP`) with SEM error bars from curve fitting.
- **Fig. 3E:** Left out — time-lapse TIRF micrographs showing single-filament life history.
- **Fig. 3F (READ — 2 readings, 34 pts):** `fig3f_binding_freq` (17 pts, 4–20 min, left y-axis 0–25 molecules) and `fig3f_p_adp_pi` (17 pts, 4–20 min, right y-axis 0–0.7, red dashed curve). Calibrated x-axis directly from the 9 printed numeric tick labels `4, 6, ..., 20`. **Key finding (`q_fig3f_adp_pi_curve_time_shift`):** the plotted red dashed `p(ADP-Pi)` curve sits consistently `9–11%` above the Melki et al. (1996) formula (`cmp_fig3f_p_adp_pi_curve` = `not reproduced`), corresponding to a `+0.65 min` (`+39 s`) horizontal shift relative to the printed `4–20 min` labels (`P_ADP-Pi(3.0 min) = 0.674` sits at the left edge `x ~ 3.65 min`).
- **Fig. 3G (READ — 2 readings, 20 pts):** `fig3g_high_adppi` (10 pts, 0.2–5 s, filled black circles, `p(ADP-Pi) > 0.5`) and `fig3g_low_adppi` (10 pts, 0.2–5 s, open white circles, `p(ADP-Pi) < 0.5`). Empirical 1-CF markers compared against single-exponential fits (`cmp_fig3g_high_adppi_curve`, `cmp_fig3g_low_adppi_caption`, `cmp_fig3g_low_adppi_panel`): `p(ADP-Pi) < 0.5` displays a small long-lived tail above `exp(-t / 1.0 s)` for `t = 1.5–5 s` (`q_fig3g_low_adppi_long_tail`).
- **Fig. 3H:** Left out — TIRF micrographs of Cy3-VASP on an aging (3–17 min) actin filament.
- **Fig. 3I (READ — 1 reading, 15 pts; fluorescence line scan left out):** `fig3i_p_adp_pi` (15 pts, 3–17 min, red dashed `p(ADP-Pi)` curve on right y-axis 0–0.7), accurately reproduced by `actin_nucleotide_state_probability` (`cmp_fig3i_p_adp_pi_curve` = `reproduced`). The noisy single-filament Cy3-VASP fluorescence intensity line scan (left y-axis, arbitrary units) was left out.

### Figure 4 (PDF p. 7 / printed p. 577)
- **Fig. 4A, 4B:** Left out — TIRF fluorescence micrographs of barbed-end elongation and Cy3-VASP tracking.
- **Fig. 4C, 4D:** Left out — TIRF kymographs (and magnified region) of single processive Cy3-VASP tetramers.
- **Fig. 4E & inset (READ — 2 readings, 52 pts):** `fig4e_be_hist` (40 histogram bars, 0.25–10.0 s) and `fig4e_inset_1cf_fit` (12 pts, 0.5–15 s, red single-exponential fit curve). Calibrated inset x-axis from tick/label positions (`x = 20` sits slightly left of right border). `cmp_fig4e_inset_1cf_curve` = `reproduced` (`tau_1 = 1.45 s`).
- **Fig. 4F:** Left out — TIRF micrograph time series of a single processive Cy3-VASP tetramer.
- **Fig. 4G (READ — 2 readings, 20 pts + SD bars):** `fig4g_measured` (9 pts, 0–100 nM VASP, black circles + SD bars) and `fig4g_mm_fit` (11 pts, 2–120 nM VASP, red Michaelis-Menten curve). Both `cmp_fig4g_mm_fit_curve` and `cmp_fig4g_measured_rates` = `reproduced`, settling `q_fig4g_unstated_vmax` (`delta_v_max = 22.9 sub/s`, `V_max = 32.9 sub/s`).
- **Fig. 4H (READ — 2 readings, 12 pts + SD bars):** `fig4h_no_vasp` (6 dark-gray bars, 0.25–2 uM actin: `1.11, 2.89, 6.19, 9.77, 12.58, 16.64 sub/s`) and `fig4h_50nM_vasp` (6 overlapping light-gray bars starting from 0: `8.53, 13.34, 20.69, 30.52, 37.39, 49.76 sub/s`), quantifying `q_fig4h_unstated_fit_intercept` (OLS slopes `8.94` and `23.5 uM^-1 s^-1`, intercepts `-0.74` and `+3.6 sub/s`).
- **Fig. 4I (READ — 2 readings, 10 pts + SEM/unstated bars):** `fig4i_dwell_time` (5 open circles + SEM bars at `0.5, 0.75, 1.0, 1.5, 2.0 uM` actin: `1.57, 1.39, 1.44, 1.28, 1.25 s`, a 20% decrease matching text `"19%"`) and `fig4i_monomers_added` (5 filled circles + unstated bars: `21.0, 32.0, 43.0, 47.0, 62.0` monomers, matching printed `42` at 1 uM and `>60` at 2 uM).

### Figure 5 (PDF p. 8 / printed p. 578)
- **Fig. 5A (READ — 9 readings across 3 subpanels, 55 pts; bottom residuals subpanels left out):** `fig5a_30uM_{green,blue,red}` (20 pts), `fig5a_60uM_{green,blue,red}` (18 pts), `fig5a_90uM_{green,blue,red}` (17 pts). Traced the three SEDPHAT best-fit absorbance-at-550-nm curves vs radial position (`5.8–7.0 cm`) in each `[hProI]` panel. Small bottom fit-residual strips were left out.
- **Fig. 5B (READ — 3 readings, 24 pts + SD bars):** `fig5b_1uM_pa`, `fig5b_1p5uM_pa`, `fig5b_2uM_pa` (8 construct conditions × 3 `[profilin-actin]` concentrations, white/gray/charcoal bars + SD bars). Key rates for condition 3 (`hProI-Actin-VASP`): `16.20, 21.8, 28.25 sub/s`; condition 4 (`hProI^H133S-Actin-VASP`): `6.30, 8.9, 11.07 sub/s`.
- **Fig. 5C (READ — 1 reading, 3 pts + SD bars):** `fig5c_100mM_kcl` (`0, 25, 50 nM` VASP at 2 uM profilin-actin, 100 mM KCl: `9.50, 12.54, 17.96 sub/s`, `1.89x` at 50 nM matching Fig. S2D(i)'s `"1.9 times faster"`).
- **Fig. 5D (READ — 2 readings, 4 pts + SEM bars):** `fig5d_hpro1` (`1.56, 1.32 s` at `1, 2 uM`) and `fig5d_hpro1_H133S` (`1.55, 1.34 s` at `1, 2 uM`).
- **Fig. 5E (READ — 2 readings, 4 pts + unstated bars):** `fig5e_hpro1` (`25.1, 37.1` monomers at `1, 2 uM`) and `fig5e_hpro1_H133S` (`10.1, 15.1` monomers at `1, 2 uM`). Quantitatively confirms `rate(5B) * tau_1(5D)` (`25.2, 37.2` and `9.7, 14.8`) and settles `q_fig5e_monomers_text_vs_plot` (p. 578 text `"~30 monomers within 1.45 s"` contradicts Fig. 5D/5E).

### Figure 6 (PDF p. 9 / printed p. 579)
- **Fig. 6A:** Left out — TIRF fluorescence micrographs of filament length in 0–20 nM MmCP ± 50 nM VASP.
- **Fig. 6B (READ — 3 readings, 22 pts):** `fig6b_no_vasp` (10 filled black circles, `0–20 nM` CP: `2873 -> 21.2` subunits), `fig6b_50nM_vasp` (9 open white circles, `0.125–20 nM` CP: `6942 -> 532` subunits), `fig6b_no_vasp_x2p4` (3 pts on dashed `-VASP x 2.4` line at `1, 2, 4 nM`: `1505, 934, 439` subunits).
- **Fig. 6C (READ — 4 readings, 32 pts):** `fig6c_4nM_no_vasp` (6 pts), `fig6c_2nM_no_vasp` (7 pts), `fig6c_4nM_50nM_vasp` (9 pts), `fig6c_2nM_50nM_vasp` (10 pts). Traced all four pure-black fit curves using the monotonic-increasing curve tracker (`trace_monotonic_capping_curve`, immune to downward-sloping callout arrows). Compared against `capping_fraction_exponential_corrected` (`1 - exp(-K_obs t)`): `2 nM -VASP` = `reproduced`, `4 nM -VASP` and `2 nM +VASP` = `close`, `4 nM +VASP` matches within `0.005–0.015`, quantitatively confirming `q_capping_eq_misprint`.
- **Fig. 6D (READ — 2 readings, 10 pts):** `fig6d_no_vasp` (5 filled circles, `0.25–4 nM`: `0.0035, 0.0090, 0.0136, 0.0235, 0.0489 s^-1`) and `fig6d_50nM_vasp` (5 open circles, `2–20 nM`: `0.0062, 0.0099, 0.0146, 0.0275, 0.0343 s^-1`). Compared against origin-constrained linear fits (`12` and `1.9 uM^-1 s^-1`) in `cmp_fig6d_no_vasp_kobs` and `cmp_fig6d_50nM_vasp_kobs`, quantifying the low-`[CP]` excess and 20 nM sublinear bend noted in `q_fig6d_sublinear_kobs_plus_vasp`.

### Figure 7 (PDF p. 10 / printed p. 580)
- **Fig. 7A, 7B:** Left out — mechanistic cartoon schematics of processive VASP barbed-end elongation and capping-protein antagonism (no numeric axes).

---

## Supplemental Information Figures (Figs. S1–S5)

### Figure S1 (PDF p. 15 / printed p. S1)
- **Fig. S1A:** Left out — cartoon schematic of VASP domain structure and GAB/FAB point mutations.
- **Fig. S1B (READ — 1 reading, 5 pts + horizontal SD bars):** `figS1b_mutants` (5 horizontal gray bars at 1 uM actin: `1uM Actin = 9.96`, `+VASP = 29.66`, `LIL-AAA = 11.83`, `RK-EE = 9.54`, `RRRK-4E = 11.67 sub/s`). Top x-axis ticks calibrated to `< 0.5 px` residual; both bar ends and horizontal SD whiskers verified on overlay.
- **Fig. S1C–S1H:** Left out — TIRF fluorescence micrographs and kymographs of wild-type and mutant Cy3-VASP binding to F-actin.

### Figure S2 (PDF p. 16 / printed p. S2)
- **Fig. S2A–S2C:** Left out — TIRF micrographs and kymographs of salt-dependent Cy3-VASP filament binding/bundling (1–100 mM KCl).
- **Fig. S2D(i) (READ — 3 readings, 15 pts):** `figS2d_i_no_vasp` (6 filled circles, `0–12 nM` CP: `2851 -> 24.8` subunits), `figS2d_i_50nM_vasp` (6 open circles, `0–12 nM` CP: `5361 -> 411` subunits; 0 nM ratio `1.88x`), `figS2d_i_no_vasp_x1p9` (3 pts on dashed `-VASP x 1.9` line at `1, 2, 4 nM`: `2047, 976, 421` subunits).
- **Fig. S2D(ii) (READ — 2 readings, 7 pts):** `figS2d_ii_no_vasp` (3 filled circles, `1, 2, 4 nM`: `0.0075, 0.0185, 0.0428 s^-1`) and `figS2d_ii_50nM_vasp` (4 open circles, `1, 2, 4, 8 nM`: `0.0040, 0.0096, 0.0133, 0.0295 s^-1`), compared against `10` and `3.7 uM^-1 s^-1` in `cmp_figS2d_ii_no_vasp_kobs` and `cmp_figS2d_ii_50nM_vasp_kobs`.
- **Fig. S2D(iii) (READ — 4 readings, 27 pts):** `figS2d_iii_{8,4,2,1}nM_no_vasp` (`6, 6, 7, 8` pts traced via `trace_monotonic_capping_curve`). Effective rates (`0.116, 0.045, 0.0191, 0.0077 s^-1`) closely match Fig. S2D(ii) `-VASP` markers (`0.0428, 0.0185, 0.0075 s^-1` at `4, 2, 1 nM`).
- **Fig. S2D(iiii) (READ — 4 readings, 26 pts):** `figS2d_iiii_{8,4,2,1}nM_vasp` (`6, 6, 7, 7` pts traced via `trace_monotonic_capping_curve`). Effective rates (`0.030, 0.0136, 0.0098, 0.0040 s^-1`) closely match Fig. S2D(ii) `+50 nM VASP` markers (`0.0295, 0.0133, 0.0096, 0.0040 s^-1`).

### Figure S3 (PDF p. 17 / printed p. S3)
- **Fig. S3A–S3C, S3F, S3G:** Left out — TIRF micrographs and kymographs of C-terminal VASP constructs, dual-color F-actin barbed-end assays, and burst-phase filament aging experiments.
- **Fig. S3D (READ — 1 reading, 4 pts + SD bars):** `figS3d_constructs` (4 vertical bars at 1.5 uM actin: `no VASP = 12.57`, `LIL-3A/RRRK-4A = 12.57`, `VASP 280-380aa = 12.25`, `GST-VASP 280-380aa = 12.76 sub/s` with SD error bars).
- **Fig. S3E (READ — 3 readings, 31 pts):** `figS3e_atp` (8 pts, 10–200 s, red solid), `figS3e_adp_pi` (12 pts, 20–1800 s, green dashed), `figS3e_adp` (11 pts, 50–2200 s, blue dotted). All three curves reproduce `actin_nucleotide_state_probability` (`k1 = 0.035 s^-1`, `k2 = 0.0026 s^-1`; `cmp_figS3e_atp_curve`, `cmp_figS3e_adp_pi_curve`, `cmp_figS3e_adp_curve` all = `reproduced`, including `P_ADP(296 s) = 0.4981` vs `0.4996`).

### Figure S4 (PDF p. 18 / printed p. S4)
- **Fig. S4A, S4B:** Left out — TIRF fluorescence micrographs of parallel and anti-parallel F-actin bundling with Cy3-VASP.
- **Fig. S4C:** Left out — single representative fluorescence line-scan intensity profile (arbitrary units, single example bundle).
- **Fig. S4D (READ — 1 reading, 3 pts + unstated bars):** `figS4d_bundling` (`single filament = 1.864`, `parallel bundle = 5.094`, `anti-parallel bundle = 5.244` fold increase Cy3-VASP signal; error bar type unstated). Quantifies `q_figs4d_fold_increase_definition` (dividing 2-filament bundle bars by 2 gives `+37%` and `+41%` per filament relative to `1.864`).

### Figure S5 (PDF pp. 19–23 / printed Fig. S5 pp. 1–5; rasterized images)
- **Fig. S5A (p. 22):** Left out — cartoon schematic of 4-site cooperative actin binding equilibria to tetrameric VASP (no numeric axes).
- **Fig. S5B (p. 23):** Left out — 3D/2D shaded MSE surface over `(K_a, beta)` (optimum `K_a = 1.11, beta = 0.2125` already recorded as scalar data).
- **Fig. S5C (p. 23; READ — 6 readings across 2 subpanels, 60 pts):** `figS5c_2uM_{5000,7000,10000}rpm` (30 pts across `X = 0.30–1.68 cm^2`) and `figS5c_25uM_{5000,7000,10000}rpm` (30 pts across `X = 0.22–1.45 cm^2`). Calibrated open L-shaped axes from inward spine ticks (`0..1.8 cm^2`, `0..12 e`-units). Quantifies the exact radial windows (`X = 0.26–1.72 cm^2`, `delta_X = 1.46 cm^2` at 2 uM actin; `X = 0.18–1.53 cm^2`, `delta_X = 1.35 cm^2` at 25 uM actin) and confirms integrated areas of `~4.4–4.8 e*cm^2` rather than unit area `1.0` (`q_auc_ka_to_kd_normalization_gap` and `q_audit_auc_25uM_x_bounds`).
