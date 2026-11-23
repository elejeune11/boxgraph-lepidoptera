# Figure Readings — `hannezo2015`

Calibrated `mode: code` digitization of all quantitative 2D plots in Hannezo, Dong, Recho, Joanny & Hayashi, *Proc. Natl. Acad. Sci. USA* **112**(28), 8620–8625 (2015) (`corpus/actin_bundling/2015_Hannezo.pdf`, pages 1–28: pp. 1–6 main article and Methods, pp. 7–25 SI Appendix Figs. S1–S16, pp. 26–28 Movie S1–S3 captions). All 31 calibrated readings were rendered at **900–1200 dpi** (plot boxes ~600–5200 px), axis frames and ticks were detected directly from the rendered pixels via `figures/_calib.py`, and every overlay PNG in `figwork/` was visually inspected after `boxgraph figure check` to verify that all magenta markers and error-bar caps lie cleanly on their intended curves, bars, or markers.

---

## Main Article (PDF pp. 1–6; PNAS 112, 8620–8625)

- **Page 1 (p. 8620, no figures):** Title, abstract, significance statement, and introduction text only; no figure panels.

- **Fig. 1A (p. 2 / p. 8621) — Partially read (`fig1a_actin_intensity`):**
  - *Read (center-top graph, `fig1a_actin_intensity`):* Relative cortical actin fluorescence intensity (%) vs. time (s) during stage-15 ring formation, 15 error-barred points across $t = 100\text{–}4250\text{ s}$ (curve ends at $t = 4362\text{ s} = 72.7\text{ min}$, matching the `72.6 min` micrograph).
  - *Calibration (1200 dpi):* X-axis from left/right axis spines ($t = 0\text{–}4500\text{ s}$, residual `0.00 px`); y-axis from 5 horizontal gridlines spanning `0–170%` (residual `3.03 px` on rasterized gridlines).
  - *Key measurement & settlement (`q_actin_increase_60_vs_70`, `q_digitized_figure_resolutions`):* Intensity rises from `104.7%` ($100\text{ s}$) through `113.9–114.9%` ($1000\text{–}1600\text{ s}$, matching SI p. 15 `"114"`) to a plateau of `153.7–156.7%` ($3700\text{–}4250\text{ s}$), confirming Main p. 1 (`"+60%"`) and SI p. 15 (`"156"`) rather than the Fig. 1A caption's `"+70%"`.
  - *Left out (left time-lapse micrographs and bottom-right stage-16 phalloidin enlargement):* Confocal microscopy images with no calibrated data axes.

- **Fig. 1B (p. 2 / p. 8621) — Left out:**
  - *Reason:* Confocal micrographs of actin (phalloidin) and myosin II heavy chain (Zip) at stage 16; no quantitative axes.

- **Fig. 1C (p. 2 / p. 8621) — Partially read (`fig1c_fourier_spectrum`):**
  - *Read (bottom-right Fourier spectrum, `fig1c_fourier_spectrum`):* Fourier amplitude vs. spatial period ($\mu\text{m}$) of early third-instar larval wild-type tracheal actin rings, 13 points across $\lambda = 0.4\text{–}2.7\ \mu\text{m}$.
  - *Calibration (1200 dpi):* Borders + interior ticks on x ($\lambda = 0\text{–}3\ \mu\text{m}$, residual `0.21 px`) and y ($0\text{–}0.10$, residual `0.23 px`).
  - *Key measurement:* Sharp peak at `1.25 um` (amplitude `0.0949`), consistent with printed `1.2 ± 0.1 um` (p. 1) and `1.2 um` (Fig. 1C caption).
  - *Left out (top micrograph and bottom-left raw intensity line scan):* Micrograph has no data axes; the raw intensity trace is a single noisy 1D fluorescence trace in arbitrary units whose sole quantitative claim is the period already captured in the Fourier spectrum.

- **Fig. 1D (p. 2 / p. 8621) — Left out:**
  - *Reason:* Conceptual schematic of the contractile actomyosin gel in frictional contact with an ECM; no quantitative axes.

- **Fig. 2A (p. 3 / p. 8622) — Partially read (`fig2a_frap_recovery`, `fig2a_frap_exp_fit`, `fig2a_frap_variance`, `fig2a_frap_variance_model`):**
  - *Read (bottom-middle recovery panel, `fig2a_frap_recovery`, `fig2a_frap_exp_fit`):* All 16 red `+` experimental markers (3 pre-bleach points at $t = 4.8\text{–}14.4\text{ s}$, intensity `0.666–0.785`; 13 post-bleach points at $t = 19.2\text{–}96.0\text{ s}$, rising from `0.020` to `0.544`) and 9 points on the green dashed single-exponential fit curve (plateau `0.561`, implying ~20% immobile fraction relative to pre-bleach baseline and single-trace recovery time $\approx 8.5\text{ s}$ vs. $n = 10$ mean $\tau = 15 \pm 7\text{ s}$).
  - *Read (bottom-right variance panel, `fig2a_frap_variance`, `fig2a_frap_variance_model`):* All 17 post-bleach error-barred experimental spatial variances $M_2(t)/M_0(t)$ ($t = 19.6\text{–}98.0\text{ s}$, `1.219` to `5.741 um^2`) and 10 points along the green dashed linear diffusion fit ($t = 10\text{–}95\text{ s}$, slope `0.0613 um^2/s` $\Rightarrow D = 0.0307\ \mu\text{m}^2/\text{s}$, intercept `-0.090 um^2`, inside the $n = 10$ population mean $D = 0.04 \pm 0.03\ \mu\text{m}^2/\text{s}$).
  - *Calibration (1200 dpi):* Bottom-middle frame + ticks (`0–100 s`, `0–0.8`), residuals `(0.37, 0.45) px`; bottom-right frame + ticks (`0–100 s`, `0–8 um^2`), residuals `(0.00, 0.19) px`.
  - *Comparisons:* Evaluated `frap_moment_recovery_eq14` using the printed population mean $D = 0.04\ \mu\text{m}^2/\text{s}$ and `var0 = 0` against both the representative trace's linear fit (`cmp_fig2a_frap_var_model`, **not reproduced**, 10/10, `max|diff| = 1.87 um^2`) and its 17 error-barred points (`cmp_fig2a_frap_var_fit`, `fit_quality`, **not reproduced**, 15/17 not reproduced, 2/17 close, `max|diff| = 2.10 um^2`), reflecting the difference between the individual representative cell ($D = 0.0307\ \mu\text{m}^2/\text{s}$) and the $n = 10$ cohort mean ($0.04 \pm 0.03\ \mu\text{m}^2/\text{s}$).
  - *Left out (top FRAP time-lapse micrographs and bottom-left spatial bleach profile):* Top images are confocal micrographs; bottom-left shows overlapping unlabelled raw fluorescence bleach profiles across multiple timepoints with arbitrary-unit intensity and unstated timestamps.

- **Fig. 2B (p. 3 / p. 8622) — Read (`fig2b_instability_threshold`):**
  - *Curve:* Red theoretical instability threshold $\chi_c(\phi) = (1 + \sqrt{\phi})^2$, 11 points across $\phi = 0.5\text{–}9.5$.
  - *Calibration (1200 dpi):* Frame borders + ticks ($\phi = 0\text{–}10$, $\chi = 0\text{–}20$), residuals `(0.33, 0.00) px`.
  - *Verdict:* Swept on `dimensionless_critical_contractility_eq10`: `cmp_fig2b_instability_threshold` **reproduced** (11/11 points, `max|diff| = 0.0175`; e.g., $\chi_c(1) = 3.983$ vs. exact `4.000`, $\chi_c(9) = 15.99$ vs. exact `16.00`).

- **Fig. 2C (p. 3 / p. 8622) — Left out:**
  - *Reason:* Confocal micrographs comparing control and 250 $\mu$M Y-27632-treated stage-16 embryos; no quantitative axes.

- **Fig. 2D (p. 3 / p. 8622) — Partially read (`fig2d_wavelength`):**
  - *Read (right graph, `fig2d_wavelength`):* Red dimensionless critical wavelength curve $\lambda_c/l_h = 2\pi\phi^{-1/4}$, 13 points across $\phi = 0.2\text{–}9.5$.
  - *Calibration (1200 dpi) & finding (`q_fig2d_shifted_interior_ticks`):* Frame borders ($\phi = 0, 10$) + tick $\phi = 2$ and y-ticks ($\lambda_c/l_h = 0\text{–}10$), residuals `(0.18, 0.39) px`. The interior x-axis tick strokes at $\phi = 4, 6, 8$ are displaced ~26.5 px (~0.21 $\phi$ units) left of their numeral labels due to a vector-graphics grouping shift, whereas borders `0, 10`, tick `2`, all digit label centers, and the plotted curve align on the true linear axis.
  - *Verdict:* Swept on `critical_wavelength_dimensionless`: `cmp_fig2d_wavelength` **reproduced** (13/13 points, `max|diff| = 0.185` on the steep $\phi = 0.2$ wing; $\phi = 1$ reads `6.282` vs. $2\pi = 6.283$).
  - *Left out (left three 1D concentration profiles):* Qualitative illustrative actin profiles with unstated axis values and no stated $\chi, \phi, \alpha$ values.

- **Fig. 2E (p. 3 / p. 8622) — Left out:**
  - *Reason:* 2D simulation spatial intensity snapshots (isotropic labyrinth vs. 20% orthoradial friction anisotropy rings) without numerical axes or colorbars (governing 2D equations also unstated; see `q_unstated_2d_anisotropic_model`).

- **Fig. 3A (p. 4 / p. 8623) — Left out:**
  - *Reason:* Confocal micrographs of wild-type, *Src42A*, and *kkv* stage-16 tracheal actin patterns; no quantitative axes.

- **Fig. 3B (p. 4 / p. 8623) — Read (`fig3b_dist_wt`, `fig3b_dist_src42`, `fig3b_bar_ring_spacing`):**
  - *Left panel (`fig3b_dist_wt`, `fig3b_dist_src42`):* Measured embryonic stage-16 ring spacing probability distributions: 13 unoccluded red WT circles ($\lambda = 0.15\text{–}1.35\ \mu\text{m}$; endpoints `0.05, 1.45` occluded at $y=0$ by green Src42 circles; distribution sum `0.997`, mean `0.494 um`, SD `0.081 um`) and all 15 green Src42A circles ($\lambda = 0.05\text{–}1.45\ \mu\text{m}$; sum `0.999`, mean `0.807 um`, SD `0.180 um`, `+63.4%` vs. WT).
  - *Right panel (`fig3b_bar_ring_spacing`):* Bar chart with broken y-axis (`0–1.2 um` lower segment, `2.0–8.0 um` upper segment) and error bars (lower halves hidden inside bars and mirrored symmetrically): `1 = Wild-type` (`0.499 ± 0.068 um`), `2 = Src42A mutant` (`0.797 ± 0.171 um`, **`+59.8%`** vs. WT), `3 = kkv mutant` (`4.575 ± 2.018 um`, **`9.18-fold`** vs. WT).
  - *Calibration (1200 dpi):* Left panel frame + ticks, residuals `(0.36, 0.29) px`; right panel lower and upper y-axis segments calibrated independently, residuals `0.33 px` and `0.30 px`.
  - *Settlement (`q_src42_and_kkv1_rounded_ratios`, `q_digitized_figure_resolutions`):* Unrounded Fig. 3B bar heights give `+59.8%` (and histograms `+63.4%`) for Src42A and `9.18-fold` for *kkv*, proving that the caption/text's `"+70%"` is an overestimate even before 1-sig-fig rounding of the means.

- **Fig. 3C (p. 4 / p. 8623) — Partially read (`fig3c_kkv_diameter_vs_actin`, `fig3c_kkv_diameter_trend`):**
  - *Read (top-right plot, `fig3c_kkv_diameter_vs_actin`, `fig3c_kkv_diameter_trend`):* All 31 red `+` scatter points showing anticorrelation between local tracheal tube diameter ($\mu\text{m}$) and local actin fluorescence concentration in *kkv1* mutants (`29.3–210.6` a.u., `5.59` down to `1.76 um`), plus 9 points on the smooth blue empirical trend curve (`5.51 um` at actin `30` to `1.85 um` at actin `205`). No theoretical equation is stated for the trend curve.
  - *Calibration (1200 dpi):* Frame borders + ticks (`0–250` a.u., `0–6 um`), residuals `(0.00, 0.47) px`.
  - *Left out (left and bottom-right panels):* Confocal micrographs of *kkv* constrictions (phalloidin and cadherin/actin staining); no quantitative axes.

- **Fig. 3D (p. 4 / p. 8623) — Left out:**
  - *Reason:* Confocal micrographs of *kkv* mutants with and without Y-27632 treatment; no quantitative axes.

- **Fig. 4A, A′, B, D, D′, E, E′, G, G′ (p. 5 / p. 8624) — Left out:**
  - *Reason:* Confocal microscopy images (live moesin-GFP imaging in control and *kkv*, 3-timepoint overlay in *kkv*, heterozygous and homozygous $\gamma$-COP mutants, and wild-type main trunk / side branch phalloidin staining); no quantitative axes.

- **Fig. 4C (p. 5 / p. 8624) — Left out:**
  - *Reason:* Qualitative 2D kymograph heatmaps comparing stationary (high friction) and chaotic (low friction) regimes with unlabelled space/time axes and no colorbar scale (quantitative nonlinear-turnover chaos phase diagram is digitized in SI Fig. S5).

- **Fig. 4F (p. 5 / p. 8624) — Read (`fig4f_cop_actin_vs_diameter`):**
  - *Points:* All 3 genotype means with both vertical and horizontal error bars (normalized cortical actin intensity vs. tube diameter in $\mu\text{m}$): `Cop-/-` (`2.511 ± 0.565 um`, `0.340 ± 0.093`; `-47.4%` diameter, `-68.9%` actin vs. WT), `Cop-/+` (`4.012 ± 0.800 um`, `0.674 ± 0.167`; `-16.0%` diameter, `-38.3%` actin vs. WT), and `Wild-type` (`4.775 ± 0.250 um`, `1.092 ± 0.144`). Refines Main p. 5 rounded percentages (`-20% / -40%` for `Cop-/+`, `-45% / -65%` for `Cop-/-`). At 1200 dpi the `Cop-/-` legend text is faint gray rather than clipped.
  - *Calibration (1200 dpi):* Frame borders + ticks (`0–6 um`, `0–1.5`), residuals `(0.20, 0.15) px`.

- **Fig. 4H (p. 5 / p. 8624) — Left out:**
  - *Reason:* Conceptual summary diagram of the proposed mechanochemical ring-formation model; no quantitative axes.

- **Page 6 (p. 8625, no figures):** Materials and Methods, Acknowledgments, and References only; no figure panels.

---

## SI Appendix (PDF pp. 7–25; SI pp. 1–19)

- **Fig. S1 (p. 7 / SI-1, panels A, A′, A″, B, B′, B″) — Left out:**
  - *Reason:* Confocal micrographs comparing taenidial folds and Lifeact-mEGFP actin rings; no quantitative axes.

- **Fig. S2 (p. 7 / SI-1) — Partially read (`figs2d_moesin_rnai_spacing`):**
  - *Read (panel D, `figs2d_moesin_rnai_spacing`):* Bar chart of late third-instar larval taenidial fold spacing ($\mu\text{m}$) with error bars: `1 = Wild-type` (`1.486 ± 0.304 um` vs. Main p. 4 `1.45 ± 0.3 um`) and `2 = Moesin RNAi` (`1.825 ± 0.517 um` vs. Main p. 4 `1.8 ± 0.5 um`; `+22.8%` increase vs. printed `25%`). (Note: Main text p. 4 swaps citations to Figs. S1 and S2; this bar chart is Fig. S2D.)
  - *Calibration (1200 dpi):* Left y-spine + 6 ticks (`0.0–2.5 um`), residual `0.00 px`.
  - *Left out (panels A, B, C):* Confocal micrographs of *moesin* RNAi embryonic and larval phenotypes; no quantitative axes.

- **Fig. S3 (p. 8 / SI-2) — Left out:**
  - *Reason:* Confocal micrographs of myosin localization in *kkv* mutants; no quantitative axes.

- **Fig. S4 (p. 8 / SI-2) — Left out:**
  - *Reason:* Confocal micrographs illustrating ImageJ ROI boxes used for actin/chitin fluorescence rescaling; no quantitative axes.

- **Pages 9–11 (SI-3–SI-5, no figures):** SI Sections 2.1–2.2 model derivations and equations only; no figure panels.

- **Fig. S5 (p. 12 / SI-6) — Read (`figs5_linear_stability`):**
  - *Curve:* Bright-red linear stability boundary $(\rho_0 \chi)_c = (1 + \sqrt{D\xi})^2$ ($D = \tau = \eta = \rho_0 = 1$, $l = 200$, $\varepsilon = 0.1$) across 11 points on the log10 friction axis $\xi \in [2 \times 10^{-5}, 0.85]$. (The discrete grid-cell phase shading for stationary vs. chaotic/transitional solutions is categorical and depends on unstated initial perturbation and $T = 5000\text{–}10^6$ classification.)
  - *Calibration (900 dpi):* Outer black box borders (`xi = 10^-5..10^0` log10, `chi*rho_0 = 0..10` linear), residuals `(0.00, 0.00) px`.
  - *Verdict:* Swept on `chaos_linear_stability_eq7`: `cmp_figs5_linear_stability` **reproduced** (11/11 points, `max|diff| = 0.0466`; e.g., $\xi = 0.1$ reads `1.734` vs. exact `1.733`, $\xi = 0.6$ reads `3.179` vs. exact `3.149`).

- **Fig. S6 (p. 12 / SI-6) — Left out:**
  - *Reason:* 2D space-time chaotic density heatmap ($\chi\rho_0 = 4, \xi = 0.001, l = 200, \varepsilon = 0.1$) from a single random initial condition; qualitative visual illustration of peak merging/emergence.

- **Fig. S7 (p. 13 / SI-7) — Left out:**
  - *Reason:* 20 densely overlapping multicolored chaotic peak-distance trajectories $\Delta_i(t)$ from a single random initial condition; individual traces are non-reproducible realization-specific paths and cross repeatedly.

- **Fig. S8 (p. 14 / SI-8) — Left out:**
  - *Reason:* Single stochastic chaotic realization of discrete peak-count trajectory $\Delta(t)$ (note caption typo prints `α = 0.1` for `ε = 0.1`); qualitative chaos illustration.

- **Page 15 (SI-9, no figures):** SI Sections 2.3–2.5 equations and zero-friction derivation only; no figure panels.

- **Fig. S9 (p. 16 / SI-10) — Read (`figs9_bifurcation_points`, `figs9_branch1`, `figs9_branch2`):**
  - *Curves:* Zero-friction ($\xi = 0$) bifurcation diagram: all 4 pitchfork bifurcation onset abscissae $\chi_i^{(0)}/\phi$ for modes $i = 1\text{–}4$ (`figs9_bifurcation_points`), plus 11 horizontal samples each along the solid linearly stable $i = 1$ branch (`figs9_branch1`) and dashed $i = 2$ branch (`figs9_branch2`) across vertical-axis values `1.2–6.5` (skipping `y = 4.0` where black dots overlap the curves). (The four small density-profile insets attached to the four black dots at $y = 4$ are illustrative profile sketches, all showing peak height `max(rho) = 4`.)
  - *Calibration (900 dpi):* Frame + interior ticks ($\chi/\phi = 0.5\text{–}4.0$, vertical axis `0–7`), residuals `(1.92, 1.68) px`.
  - *Verdicts & TWO major discoveries (`q_figs9_yaxis_maxrho_not_l2`, settles `q_si_section25_xi0_curvature_and_figs9`):*
    1. **Captioned $\alpha^2\phi = 10$ refuted; true plotted parameter is $\alpha^2\phi = 100$:** Digitized onset abscissae `1.093, 1.387, 1.881, 2.572` refute `cmp_figs9_onsets_captioned_10` (**not reproduced**, 4/4, predicts `1.987..16.79`, `max|diff| = 14.22`) and match `cmp_figs9_onsets_100` (**reproduced**, 4/4, predicts `1.099, 1.395, 1.888, 2.579`, `max|diff| = 0.0076`).
    2. **Vertical axis mislabeled as $\|\rho\|_2$ — actually plots $\max(\rho)$ ($L^\infty$ norm):** Exact harmonic continuation of the SI Section 2.5 ODE at $\alpha^2\phi = 100$ refutes the captioned $L^2$ norm (`zero_friction_ode_branch_l2`: `cmp_figs9_branch1_l2` **not reproduced**, 11/11, `max|diff| = 2.11`; `cmp_figs9_branch2_l2` **not reproduced**, 11/11, `max|diff| = 5.72`) and reproduces both branches to `< 0.25%` when evaluated as $\max(\rho)$ (`zero_friction_ode_branch_maxrho`: `cmp_figs9_branch1_maxrho` **reproduced**, 11/11, `max|diff| = 0.0031`; `cmp_figs9_branch2_maxrho` **reproduced**, 11/11, `max|diff| = 0.0026`). Both `ref_figs9_branch1_modes` (`rel_chg = 0.000000`) and `ref_figs9_branch1_l2_modes` (`rel_chg = 0.000134`) settle cleanly below `0.1%`.

- **Fig. S10 (p. 17 / SI-11) — Read (`figs10_finite`, `figs10_infinite`):**
  - *Curves:* Number of patterns vs. $\phi$ at $\alpha = 25$: blue finite-domain step curve $i_c/2$ across all 6 plateaus (`figs10_finite`, 6 points at $\phi = 0.07, 0.16, 0.33, 0.70, 1.05, 1.65$) and purple infinite-domain curve $\alpha\phi^{1/4}/(2\pi)$ (`figs10_infinite`, 10 points across $\phi = 0.05\text{–}1.95$).
  - *Calibration (900 dpi):* Frame borders + ticks ($\phi = 0\text{–}2$, $n = 0\text{–}5$), residuals `(0.35, 2.40) px`.
  - *Verdicts:* `cmp_figs10_infinite` **reproduced** (10/10 points, `max|diff| = 0.0262`); `cmp_figs10_finite` **close** (5/6 points reproduced, 1/6 close at $\phi = 1.65$: exact `4.5` vs. digitized `4.487 [4.475, 4.498]`, diff `0.013`, 0.002 outside the narrow 7 px stroke band due to the 2.40 px y-tick scatter).

- **Page 18 (SI-12, no figures):** SI Section 2.6.2 weakly nonlinear perturbation derivation only; no figure panels.

- **Fig. S11 (p. 19 / SI-13) — Read (`figs11_bifurcation_points`, `figs11_primary_branch`):**
  - *Curves:* Subcritical AUTO bifurcation diagram at $\alpha = 35, \phi = 0.1$: all 7 red `B` linear bifurcation thresholds $\chi_i^{(0)}$ along $\max\rho = 1$ for modes $i = 4\text{–}10$ (`figs11_bifurcation_points`), plus 13 horizontal samples along the primary $i_c = 6$ branch (`figs11_primary_branch`: lower subcritical leg $\max\rho = 1.2\text{–}3.2$ and upper branch $\max\rho = 6.5\text{–}10.5$, noting that above $\max\rho \approx 3.8$ a secondary bifurcation branches left of the pure $i = 6$ curve so the primary branch is the 2nd curve from the left).
  - *Calibration (900 dpi):* Frame + ticks ($\chi = 1.60\text{–}2.05$, $\max\rho = 0\text{–}11$), residuals `(1.75, 1.73) px`.
  - *Verdicts:*
    - `cmp_figs11_bif_points` (`finite_domain_bifurcation_threshold`): **reproduced** (7/7 points, `max|diff| = 0.000109`).
    - `cmp_figs11_branch_exact` (`steady_state_eq8_branch_continuation`): **reproduced** (13/13 points across both lower subcritical and upper legs, `max|diff| = 0.00136`; spectral refinements `ref_figs11_branch_exact_modes_chi` `rel_chg = 0.000017` and `ref_figs11_branch_exact_modes_minrho` `rel_chg = 0.000420` both settle `< 0.1%`).
    - Weakly nonlinear quadratic approximations (`purpose: approximation`): both `cmp_figs11_branch_wnl_printed` (**not reproduced**, 8/13 within `phi = 0.1` rounding band, 5/13 not reproduced on upper leg) and `cmp_figs11_branch_wnl_corrected` (**not reproduced**, 5/13, 8/13) inevitably fail past the fold ($\max\rho \gtrsim 3$). Near onset ($\max\rho = 1.2$, exact $\chi = 1.7296$, digitized `1.730`), printed $\chi_i^{(2)} = -0.0384$ barely moves (`1.7341`, 10× too flat) whereas adjoint-Fredholm $\chi_i^{(2)} = -0.4005$ captures the onset curvature (`1.7268`, before the 2nd-harmonic shift between $\varepsilon$ and $\max\rho - 1$ grows at larger amplitude), confirming `q_si_section262_wrong_fredholm_eigenvector`.

- **Fig. S12 (p. 20 / SI-14) — Read (`figs12_bifurcation_points`, `figs12_primary_branch`):**
  - *Curves:* Supercritical AUTO bifurcation diagram at $\alpha = 15, \phi = 3$: 8 isolated red `B` linear bifurcation thresholds for modes $i = 4, 6, 7, 9, 10, 11, 12, 13$ (`figs12_bifurcation_points`; the overlapping $i = 5, 8$ glyph pair is excluded cleanly by glyph width), plus 10 horizontal samples along the primary $i_c = 6$ branch (`figs12_primary_branch`, $\max\rho = 1.3\text{–}13.0$).
  - *Calibration (900 dpi):* Frame + ticks ($\chi = 7\text{–}14$, $\max\rho = 0\text{–}14$), residuals `(1.56, 1.61) px`.
  - *Verdicts:*
    - `cmp_figs12_bif_points` (`finite_domain_bifurcation_threshold`): **reproduced** (8/8 points, `max|diff| = 0.00457`, e.g., $i = 6$ reads `7.480` vs. `7.479`, $i = 13$ reads `11.82` vs. `11.818`).
    - `cmp_figs12_branch_exact` (`steady_state_eq8_branch_continuation`, `M_modes: 60`): **reproduced** (10/10 points across $\max\rho = 1.3\text{–}13.0$, `max|diff| = 0.00837`; `ref_figs12_branch_exact_modes_chi` `rel_chg = 0.000006` and `ref_figs12_branch_exact_modes_minrho` `rel_chg = 0.000138` both settle `< 0.1%`).
    - Weakly nonlinear quadratic approximations (`purpose: approximation`): `cmp_figs12_branch_wnl_printed` (**not reproduced**, 5/10, `max|diff| = 34.1`) and `cmp_figs12_branch_wnl_corrected` (**not reproduced**, 3/10, `max|diff| = 248.5`) both diverge quadratically at large amplitude ($\max\rho \ge 3.5$), while near onset ($\max\rho = 1.3$, digitized `7.556`) printed $\chi_i^{(2)} = +0.5546$ is 6.4× too shallow (`7.504`) and adjoint-Fredholm $\chi_i^{(2)} = +3.5325$ gives `7.638` (overshooting slightly because $\max\rho - 1 > \varepsilon$ due to the positive 2nd harmonic).

- **Page 21 (SI-15, no figures):** SI Sections 2.6.4–2.6.5 parameter estimation and simulation description only; no figure panels.

- **Fig. S13 (p. 22 / SI-16) — Read (`figs13_S_vs_t`, `figs13_chi_ramp`):**
  - *Main plot (`figs13_S_vs_t`):* Ensemble-averaged order parameter $S(t)$, 12 points across $t = 100\text{–}2650\text{ s}$ (`45.4` at $100\text{ s}$, jumping near $t \approx 170\text{ s}$ to `1.038e4` at $250\text{ s}$, plateauing at `1.392e4` for $t \ge 2000\text{ s}$ — ~116× larger than continuum $\int_0^1 (\rho - 1)^2 dy \approx 120$, confirming `q_unstated_pde_ic_and_subcritical_sensitivity`). Recorded as `origin: reported` due to unstated initial noise spectrum and normalization/grid factor.
  - *Inset (`figs13_chi_ramp`):* Imposed linear Péclet-number ramp $\chi(t)$, 7 points across $t = 100\text{–}2500\text{ s}$ (`1.099` at $100\text{ s}$ to plateau `1.475` for $t \ge 2100\text{ s}$).
  - *Calibration (1200 dpi):* Main frame + ticks ($t = 0\text{–}3000\text{ s}$, $S = 0\text{–}1.6 \times 10^4$), residuals `(2.68, 2.02) px`; inset frame + ticks ($t = 0\text{–}3000\text{ s}$, $\chi = 1.0\text{–}1.6$), residuals `(0.00, 2.37) px`.
  - *Verdicts (inset):* `cmp_figs13_chi_ramp_printed` (`linear_chi_ramp_si` with printed rounded `1.07 -> 1.47` over `1875 s`) is **close** (6/7 reproduced, 1/7 close at $t = 100\text{ s}$: `1.0913` vs. `1.099 [1.092, 1.107]`); `cmp_figs13_chi_ramp_chained` (using unrounded Eq. (10) `chi_c = 1.076` and chained plateau `1.473`) is **reproduced** (7/7 points, `max|diff| = 0.00583`), showing the inset ramp used the unrounded Eq. (10) values.

- **Fig. S14 (p. 22 / SI-16) — Read (`figs14_wavelength_histogram`):**
  - *Histogram:* All 28 non-zero histogram bars of steady-state simulated ring wavelengths ($N = 80$ realizations, $\lambda = 0.529\text{–}1.491\ \mu\text{m}$), peaking at `0.608 um` (probability density `12.33 um^-1`) with a secondary period-doubling bump around `1.04–1.12 um` (histogram-weighted mean `~0.68 um`, consistent with $10\ \mu\text{m} / 14.68$ rings). Recorded as `origin: reported` because initial noise spectrum and discretization are unstated.
  - *Calibration (1200 dpi):* X-axis from frame + ticks (`0.5–1.5 um`, residual `1.83 px`); y-axis from top border + 6 interior ticks (`2–14 um^-1`, omitting bar-thickened bottom border, residual `0.23 px`).

- **Fig. S15 (p. 22 / SI-16) — Read (`figs15_n_rings_vs_chi`, `figs15_sqrt_fit`):**
  - *Points & curve:* All 38 error-barred ensemble-mean steady-state ring counts (`figs15_n_rings_vs_chi`, $\chi = 1.11\text{–}5.00$, rising from `10.49 ± 1.22` through `14.87 ± 0.87` at $\chi = 1.42$ and `15.84 ± 1.04` at $\chi = 1.53$ to `34.35 ± 1.40` at $\chi = 5.00$) and 9 points along the red square-root fit curve (`figs15_sqrt_fit`, $\chi = 1.2\text{–}4.9$, closely following $n \approx 15.5\sqrt{\chi}$; no fit equation is printed in the paper). Recorded as `origin: reported` / `fitted` due to unstated ensemble IC/discretization and unprinted fit prefactor.
  - *Calibration (1200 dpi):* X-axis from frame + ticks ($\chi = 0\text{–}6$, residual `0.30 px`); y-axis from 7 interior ticks ($n = 5\text{–}35$, omitting asterisk-thickened $n = 0$ bottom border, residual `0.23 px`).

- **Page 23 (SI-17, no figures):** SI Section 2.7 two-cell coordination equations (Eqs. 15–16) only; no figure panels.

- **Fig. S16 (p. 24 / SI-18) — Partially read (`figs16_correlation_vs_beta`):**
  - *Read (top panel, `figs16_correlation_vs_beta`):* All 13 blue markers of intercellular pattern correlation $C(\beta)$ across $\beta = 1.02 \times 10^{-5}\text{–}9.87 \times 10^{-2}$: $C = 0.0131 \pm 0.0022$ at $\beta = 1.02 \times 10^{-5}$, crossing $C = 0.5023$ at $\beta = 1.00 \times 10^{-3}$, reaching $C = 0.948$ at $\beta = 5.11 \times 10^{-3}$, and saturating at $C = 1.001$ for $\beta \ge 0.05$. Recorded as `origin: reported` because Eq. (16) omits spatial noise correlation and cell parameters.
  - *Calibration (900 dpi):* Log10 x-axis ticks (`10^-6..10^0`, residual `0.41 px`) and linear y-axis ticks (`0.0–1.2`, residual `0.41 px`).
  - *Settlement (`q_unstated_two_cell_coupling_and_C_formula`, `q_digitized_figure_resolutions`):* Because positive densities oscillating in $[0.6, 1.5]$ around $\rho_0 = 1$ (bottom panels) would yield $C_{\text{uncorr}} = 1/(1 + \mathrm{Var}(\rho)) \approx 0.90$ under the printed formula $C = \int \rho_1 \rho_2 / \int \rho_1^2$, the measured $C = 0.0131$ at $\beta = 1.02 \times 10^{-5}$ proves mean subtraction $(\rho_i - \rho_0)$ was used in code despite being omitted in print.
  - *Left out (bottom three panels):* Three single-run spatial density snapshots ($\rho_1, \rho_2$ on $x \in [0, 400]$ at $\beta = 0.005$) at unstated timepoints and unstated PDE parameters; qualitatively illustrate alignment dynamics.

- **Page 25 (SI-19, no figures):** SI Appendix References only; no figure panels.

---

## Supporting Movies (PDF pp. 26–28)

- **Movie S1 (p. 26), Movie S2 (p. 27), Movie S3 (p. 27) — Left out:**
  - *Reason:* Single-frame confocal video thumbnails and captions (stage-15 actin ring live imaging, FRAP live imaging, and *kkv1* moesin-GFP live imaging); no quantitative axes.
- **Page 28 (no figures):** List of other Supporting Information files only; no figure panels.
