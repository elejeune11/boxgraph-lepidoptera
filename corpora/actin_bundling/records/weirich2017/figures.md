# Figure Inventory and Digitization Log — `weirich2017`

Complete panel-by-panel audit of every page of `corpus/actin_bundling/2017_Weirich.pdf` (16 PDF pages: pp. 1–6 main PNAS article printed pp. 2131–2136; pp. 7–16 Supporting Information printed pp. SI-1–SI-10). All 40 quantitative figure readings across all 13 quantitative panels were implemented as 1200-dpi `mode: code` scripts under `figures/`, verified through `boxgraph figure check`, and visually inspected via their magenta overlay PNGs in `figwork/`.

## Page 1 (printed p. 2131)
- No figures (Title, Abstract, Significance, and Introduction text only).

## Page 2 (printed p. 2132) — Figure 1: Cross-linked F-actin liquid-crystal droplets and phase diagram
- **Fig. 1A**: Left out — 3 fluorescence micrographs of TMR-actin at 0, 1, and 5 mol % filamin (`10 um` scale bar); qualitative morphology images with no quantitative axes.
- **Fig. 1B**: Left out — 3 fluorescence micrographs of a tactoid (TMR-actin, Alexa-488 filamin, and merged composite; `5 um` scale bar); qualitative colocalization microscopy.
- **Fig. 1C (top strip — FRAP micrographs)**: Left out — 4 fluorescence time-lapse micrographs (`0, 5, 10, 15 min`, `2 um` scale bar) illustrating photobleaching recovery in a single tactoid.
- **Fig. 1C (bottom plot — FRAP recovery curve)**: Digitized at 1200 dpi (`cal_res = (0.70, 0.62) px`):
  - `fig1c_frap_data` (`d_fig1c_frap_data`): all 14 black filled circles (`t = -0.201` pre-bleach at `I = 0.932`, plus 13 post-bleach frames `t = 0.05..16.48 min`, `I = 0.223..0.774`, reproducing the 5.05–8.15 min sampling gap also seen in SI Fig. S1).
  - `fig1c_frap_fit` (`d_fig1c_frap_fit`): 8 samples along the black dashed exponential fit (`t = 1..15 min`, `I = 0.282..0.742`). An unconstrained 3-parameter exponential fit to the digitized dashed curve gives `I_0 = 0.224, I_inf = 0.947, tau = 712 s` (caption prints `tau_R = 880 s`, which corresponds to fixing `I_inf ~ 1.03–1.05`; logged in `q_digitized_fig3a_slope_and_fig1c_fit_tau`).
- **Fig. 1D (actin/filamin phase diagram: filamin mole fraction vs filament length)**: Digitized all 4 phase marker classes at 1200 dpi (`cal_res = (0.65, 0.60) px`, calibrating `x` on all 19 log ticks `[0.1..10.0]` because the right box spine sits slightly past `10 um`):
  - `fig1d_gas` (`d_fig1d_gas`): 10 black `+` markers (`L = 0.121..0.543 um`, `filamin = 0..1%`).
  - `fig1d_liquid_droplets` (`d_fig1d_liquid_droplets`): 17 cyan filled circles (`L = 0.121..0.546 um`, `filamin = 1..10%`).
  - `fig1d_fluid_bundles` (`d_fig1d_fluid_bundles`): 12 dark-blue open circles (`L = 0.364..1.085 um`, `filamin = 2.5..10%`).
  - `fig1d_viscoelastic_networks` (`d_fig1d_viscoelastic_networks`): 12 gray `x` markers (`L = 0.778..10.04 um`, `filamin = 0..10%`).
  - Qualitative dashed phase-boundary guide curves left out.

## Page 3 (printed p. 2133) — Figure 2: Cross-link density tunes tactoid aspect ratio and anchoring strength
- **Fig. 2A (top strip — tactoid micrographs)**: Left out — 4 fluorescence micrographs of representative tactoids at 2.5, 5, 10, and 15 mol % filamin (`5 um` scale bar).
- **Fig. 2A (bottom plot — dual-axis aspect ratio $L/r$ and anchoring strength $\omega$ vs filamin concentration)**: Digitized both series with their full vertical error bars at 1200 dpi (`cal_res = (0.48, 0.14) px`, `error_bars: unstated` because caption does not specify SD vs SEM):
  - `fig2a_aspect_ratio` (`d_fig2a_aspect_ratio`): 4 black open circles at 2.5, 5, 10, 15 mol % filamin: `L/r = 2.938 [2.691, 3.191], 2.808 [2.582, 3.041], 2.438 [2.241, 2.640], 2.037 [1.750, 2.331]`.
  - `fig2a_omega` (`d_fig2a_omega`): 4 green filled diamonds: `omega = 2.164 [1.813, 2.551], 1.977 [1.669, 2.316], 1.488 [1.257, 1.747], 1.039 [0.765, 1.359]`.
  - **Model comparison**: `cmp_fig2a_omega_from_ar` (`anchoring_strength_from_aspect_ratio_curve`, $\omega = (L/r)^2/4$) -> **`reproduced`** (4/4 points `reproduced`, `2.158, 1.971, 1.486, 1.038` vs `2.164, 1.977, 1.488, 1.039`, matching to 0.2% including the exact asymmetric error-bar propagation `(L/r ± sigma)^2 / 4`). Also notes that the 15% point sits at `omega = 1.04`, slightly below SI p. 8's printed `1.1 <= omega <= 2.2`.
- **Fig. 2B**: Left out — fluorescence micrograph, PolScope polarization retardance/orientation image, and nematic director schematic (`2 um` scale bar); no quantitative coordinate axes.
- **Fig. 2C**: Left out — qualitative conceptual shape diagram of tactoids on unnumbered $\gamma_A$ vs $\gamma_I$ axes (no numeric ticks or scale).

## Page 4 (printed p. 2134) — Figure 3: Tactoid coarsening, coalescence dynamics, and continuum model comparison
- **Fig. 3A (top strip — coarsening micrographs)**: Left out — 4 fluorescence time-lapse micrographs (`10', 25', 50', 85'`, `5 um` scale bar) of tactoid growth.
- **Fig. 3A (bottom plot — dual-axis log-log average length and aspect ratio vs time)**: Digitized all 4 series at 1200 dpi (`cal_res = (0.74, 1.32) px` left length axis, `(0.74, 0.64) px` right aspect-ratio axis; `error_bars: sd` per caption):
  - `fig3a_length_data` (`d_fig3a_length_data`): 16 black closed circles with ±1 SD error bars (`t = 10.0..85.4 min`, `L = 1.461..4.237 um`, single-dataset log-log OLS slope `0.4920`).
  - `fig3a_length_fit` (`d_fig3a_length_fit`): 8 points along the black dashed power-law fit line (`t = 12..60 min`, `L = 1.631..3.616 um`, slope `0.4936` — matching this individual dataset and close to theory `0.50`, compared to the caption's 4-dataset mean `0.47 ± 0.01`; logged in `q_digitized_fig3a_slope_and_fig1c_fit_tau`).
  - `fig3a_aspect_ratio_data` (`d_fig3a_aspect_ratio_data`): 16 gray open squares with ±1 SD bars (`L/r = 2.037..2.163`).
  - `fig3a_aspect_ratio_fit` (`d_fig3a_aspect_ratio_fit`): 8 samples along the near-constant solid gray trend line (`L/r = 2.071..2.113`).
  - **Model comparisons**: `cmp_fig3a_fit_measured_alpha` (`reproduced`, 8/8, prefactor `0.5173`), `cmp_fig3a_fit_theory_alpha` (`reproduced`, 8/8, prefactor `0.4683`, 4x smaller max residual `0.016 um` than `alpha = 0.47`), `cmp_fig3a_data_measured_alpha` (`reproduced`, 16/16), `cmp_fig3a_data_theory_alpha` (`reproduced`, 16/16).
- **Fig. 3B (top strip — coalescence micrographs)**: Left out — 4 fluorescence time-lapse micrographs (`-10 s, 0 s, 30 s, 125 s`, `5 um` scale bar) of two coalescing tactoids.
- **Fig. 3B (bottom plot — major-axis length vs time during coalescence)**: Digitized at 1200 dpi (`cal_res = (0.44, 0.55) px`):
  - `fig3b_pre_upper` (`d_fig3b_pre_upper`): 10 interior pre-coalescence blue dots of the larger tactoid (`t = -89.2..+0.8 s`, `L = 5.10..5.36 um`).
  - `fig3b_pre_lower` (`d_fig3b_pre_lower`): 10 pre-coalescence blue dots of the smaller tactoid (`t = -94.3..-4.3 s`, `L = 4.06..4.36 um`).
  - `fig3b_coalesced_data` (`d_fig3b_coalesced_data`): 48 consecutive 5-s post-coalescence blue dots (`t = 10.8..245.8 s`, `L = 9.25` decaying to `5.76..6.17 um`).
  - `fig3b_exp_fit` (`d_fig3b_exp_fit`): 12 samples along the black dashed fit curve (`t = 15..240 s`, `L = 8.789..6.011 um`; exponential fit gives `L_f = 6.005 um`, `L(10 s) = 9.198 um`, `tau = 36.29 s`, matching the Fig. 3D 10% point at `(120 nm/s, 36.4 s)`).
  - Inset grey tactoid cartoons labeled $L_i$ and $L_f$ left out.
- **Fig. 3C (main panel — rescaled coalescence collapse)**:
  - `fig3c_master_curve` (`d_fig3c_master_curve`): 12 samples along the solid black exponential decay curve at 1200 dpi (`cal_res = (0.20, 1.78) px`, `t/tau = 0.2..11.5`, `y = 0.812..-0.0025`).
  - **Model comparison**: `cmp_fig3c_master_curve` (`tactoid_exponential_master_curve`, $\exp(-t/\tau)$) -> **`reproduced`** (12/12 points `reproduced`, max error `0.0069`).
  - Overplotted multi-color scatter of 14 overlapping individual coalescence trajectories left out (their fitted parameters are digitized individually in the inset below).
- **Fig. 3C (Inset — $\tau$ vs $L_i - L_f$)**: Digitized all 3 series at 1200 dpi (`cal_res = (0.30, 2.32) px`):
  - `fig3c_inset_5pct` (`d_fig3c_inset_5pct`): all 8 red crosses (`L_i - L_f = 0.81..2.27 um`, `tau = 8.30..16.30 s`).
  - `fig3c_inset_10pct` (`d_fig3c_inset_10pct`): all 6 blue circles (`L_i - L_f = 1.86..4.37 um`, `tau = 24.53..58.51 s`).
  - `fig3c_inset_linear_fit` (`d_fig3c_inset_linear_fit`): 7 points along the dashed black linear regression line (`L_i - L_f = 1.0..4.0 um`, `tau = 5.39..46.56 s`), corresponding to $\tau = 13.72\,(L_i - L_f) - 8.34\text{ s}$.
  - **Model comparison**: `cmp_fig3c_inset_linear_fit` (`coalescence_tau_linear_scaling`, zero-intercept $\tau = (\eta/\gamma_{\text{eff}})(L_i - L_f)$ at slope `10 s/um`) -> **`not reproduced`** (5/7 `not reproduced`, 2/7 `close`), confirming the plotted regression has a substantial negative intercept (`-8.3 s`) and steeper slope (`13.7 s/um`; logged in `q_digitized_fig3c_inset_nonzero_intercept`).
- **Fig. 3D ($\tau$ vs $\gamma_I/\eta$)**: Digitized all 4 series at 1200 dpi (`cal_res = (0.30, 0.30) px`):
  - `fig3d_5pct_data` (`d_fig3d_5pct_data`): all 8 red crosses (`gamma_I/eta = 100.0..160.1 nm/s`, `tau = 8.30..15.50 s`, resolving the ~23-px vertically stacked pair at `160.0 nm/s` via template matching pursuit + XOR refinement).
  - `fig3d_10pct_data` (`d_fig3d_10pct_data`): all 6 blue circles (`gamma_I/eta = 83.0..160.0 nm/s`, `tau = 24.60..56.60 s`).
  - `fig3d_theory_w2p0` (`d_fig3d_theory_w2p0`): 9 points along the lower dashed theory curve labeled $\omega = 2.0$ (`gamma_I/eta = 25..170 nm/s`, `tau = 51.45..7.51 s`; exact hyperbola $(\gamma_I/\eta)\tau = 1278 \pm 2\text{ nm}$).
  - `fig3d_theory_w1p4` (`d_fig3d_theory_w1p4`): 8 points along the upper dashed theory curve labeled $\omega = 1.4$ (`gamma_I/eta = 60..175 nm/s`, `tau = 60.09..20.63 s`; exact hyperbola $(\gamma_I/\eta)\tau = 3607 \pm 3\text{ nm}$, a **2.822x** separation above $\omega = 2.0$).
  - **Model comparisons (settles BOTH the Rayleigh factor of 2 and `q_fig3d_omega_vs_volume_separation`)**:
    - `cmp_fig3d_w2p0_printed` (`tactoid_shape_dynamics_ode`, printed `D = (1/2) V eta [...]`): **`not reproduced`** (9/9, `614 nm` vs `1278 nm`, **2.08x too fast**).
    - `cmp_fig3d_w2p0_rayleigh2` (`tactoid_shape_dynamics_ode_rayleigh2`, `P = 2*Phi`): **`reproduced`** (9/9 `reproduced`, `1228 nm` vs `1278 nm`, within 3.9% and inside Table S1 `L_eq` uncertainty band).
    - `cmp_fig3d_w1p4_printed`: **`not reproduced`** (8/8, `1577 nm` vs `3607 nm`, **2.29x too fast**).
    - `cmp_fig3d_w1p4_rayleigh2`: **`not reproduced`** (8/8, `3153 nm` vs `3607 nm`, systematically 12.6% low at Table S1 mean `L_eq = 6.04 um`, matching closely if `L_eq ~ 6.9 um`).
    - `cmp_fig3d_w1p4_fixedV_rayleigh2`: **`not reproduced`** (8/8, holding volume fixed at the 5% droplet volume `V = 1.088 um^3`, `L_eq = 2.455 um` changes `v_cap*tau` by only +5.8% to `1299 nm`, **2.78x below** `3607 nm` — directly disproving SI p. 9's claim that the two Fig. 3D curves share the same volume; logged in `q_digitized_fig3d_rayleigh2_and_volume_settled`).

## Page 5 (printed p. 2135) — Figure 4: Capillary pearling and spontaneous contraction of fluid bundles
- **Fig. 4A**: Left out — 3 fluorescence time-lapse micrographs (`0, 10, 40 min`, `10 um` scale bar) of bundle pearling instability.
- **Fig. 4B**: Left out — 4 fluorescence time-lapse micrographs (`0, 2, 5, 10 min`, `10 um` scale bar) of spontaneous bundle contraction.
- **Fig. 4C (normalized bundle length $L(t)/L_i$ vs time)**: Digitized all 4 series at 1200 dpi (`cal_res = (0.00, 0.15) px`):
  - `fig4c_2p5pct_data` (`d_fig4c_2p5pct_data`): 9 binned positions (`t = 50..1450 s`, `L/L_i = 0.964..0.717`) through the dense red open-circle trajectory.
  - `fig4c_2p5pct_fit` (`d_fig4c_2p5pct_fit`): 9 samples along the upper black dashed continuum fit (`t = 100..1450 s`, `L/L_i = 0.937..0.730`; exponential fit yields `L_f/L_i = 0.7242`, `tau = 375.6 s`, RMS `0.00026`).
  - `fig4c_10pct_data` (`d_fig4c_10pct_data`): 9 binned positions (`t = 40..1050 s`, `L/L_i = 0.906..0.550`) through the dense blue filled-circle trajectory.
  - `fig4c_10pct_fit` (`d_fig4c_10pct_fit`): 9 samples along the lower black dashed continuum fit (`t = 100..1450 s`, `L/L_i = 0.823..0.540`; exponential fit yields `L_f/L_i = 0.5390`, `tau = 199.0 s`, RMS `0.00086`).
  - Inset grey bundle cartoons labeled $L_i$ and $L_f$ left out.
- **Fig. 4D (rescaled bundle contraction master collapse)**: Digitized binned medians and 10–90% spread envelopes for all 3 filamin concentrations at 1200 dpi (`cal_res = (0.30, 0.22) px`):
  - `fig4d_2p5pct` (`d_fig4d_2p5pct`): 7 bins (`t/tau = 0.25..3.0`, green points) -> `cmp_fig4d_2p5pct_collapse`: **`close`** (6/7 `reproduced`, 1 `close`, max diff `0.029`).
  - `fig4d_5pct` (`d_fig4d_5pct`): 9 bins (`t/tau = 0.25..4.0`, red open circles) -> `cmp_fig4d_5pct_collapse`: **`reproduced`** (9/9 `reproduced`, max diff `0.036`).
  - `fig4d_10pct` (`d_fig4d_10pct`): 7 bins (`t/tau = 0.25..3.0`, blue dots) -> `cmp_fig4d_10pct_collapse`: **`reproduced`** (7/7 `reproduced`, max diff `0.060`).
- **Fig. 4E (bundle contractile strain box-and-whisker plot vs filamin concentration)**:
  - `fig4e_contractile_strain` (`d_fig4e_contractile_strain`, 1200 dpi, `cal_res_y = 0.96 px`, `error_bars: unstated`): calibrated strictly from the 4 unambiguous regular ticks `[0.2, 0.4, 0.5, 0.6]` while excluding the misprinted `0.6` tick at the `0.3` position (which is also physically shifted ~8 px upward from the regular `407.6 px` tick grid).
  - Extracted medians (box Q1–Q3; whiskers): `2.5%`: `0.271` (`[0.233, 0.278]`; `[0.220, 0.281]`); `5.0%`: `0.350` (`[0.330, 0.423]`; `[0.281, 0.460]`); `10.0%`: `0.520` (`[0.482, 0.536]`; `[0.400, 0.581]`).

## Page 6 (printed p. 2136)
- No figures (Discussion, Materials and Methods, Acknowledgments, and References text only).

## Pages 7–8 (printed pp. SI-1–SI-2)
- No figures (SI Materials and Methods and Equilibrium Droplet Shapes derivation text only).

## Page 9 (printed p. SI-3) — SI Figure S1: Raw FRAP intensity traces
- **Fig. S1 (main plot — bleached and unbleached fluorescence intensity vs time)**: Digitized both dense dot series across `t = 0..800 s` (respecting the `330..480 s` acquisition gap) at 1200 dpi (`cal_res = (0.40, 0.47) px`):
  - `figS1_bleached` (`d_figS1_bleached`): 13 post-bleach binned positions along the pink dots (`I = 230.5..646.4`).
  - `figS1_unbleached` (`d_figS1_unbleached`): 13 post-bleach binned positions along the blue dots (`I = 1043..864..940`). Their bleached/unbleached ratio matches Fig. 1C across the entire 0–800 s window.
- **Fig. S1 (bottom-right inset micrograph)**: Left out — fluorescence micrograph of a photobleached tactoid indicating the blue (unbleached) and pink (bleached) ROIs; no plot axes.

## Page 10 (printed p. SI-4) — SI Figures S2 and S3: Continuum tactoid free energy, aspect ratio, and instability geometry
- **Fig. S2A**: Left out — geometric schematic of a homogeneous tactoid defining surface normal $\hat{N}$ and nematic director field $\hat{n}$.
- **Fig. S2B**: Left out — geometric schematic defining circular-arc tactoid parameters $L$, $r$, $R$, and $\alpha$.
- **Fig. S2C (rescaled free energy $F/(\gamma_I V^{2/3})$ vs $\alpha$)**: Digitized all 3 curves at 1200 dpi (`cal_res = (0.23, 0.20) px`):
  - `figS2c_Fhat_w0` (`d_figS2c_Fhat_w0`, blue, 9 pts `alpha = 0.3..3.0`, min `2.621` at `alpha = pi/2`).
  - `figS2c_Fhat_w0p5` (`d_figS2c_Fhat_w0p5`, orange, 9 pts `alpha = 0.3..3.0`, min `2.974` at `alpha = 1.1`).
  - `figS2c_Fhat_w2` (`d_figS2c_Fhat_w2`, green, 8 pts `alpha = 0.2..3.0`, min `3.522` at `alpha = 0.67`).
  - Inset colored droplet schematics left out.
  - **Model comparisons (settles `q_figS2C_missing_2pi_one_third_factor`)**:
    - Exact `tactoid_rescaled_free_energy` ($F/(\gamma_I V^{2/3})$): `cmp_figS2c_Fhat_w0`, `cmp_figS2c_Fhat_w0p5`, `cmp_figS2c_Fhat_w2` -> all 3 **`not reproduced`** (26/26 points systematically higher by an exact constant ratio `1.8453 ± 0.001 = (2*pi)^(1/3)`).
    - Corrected variant `tactoid_rescaled_free_energy_figS2c_norm` ($F/[(2\pi)^{1/3}\gamma_I V^{2/3}] = (2f + \omega g)/g^{2/3}$): `cmp_figS2c_Fhat_norm_w0` (9/9), `cmp_figS2c_Fhat_norm_w0p5` (9/9), `cmp_figS2c_Fhat_norm_w2` (8/8) -> all 3 **`reproduced`** (26/26 points, max error `0.001–0.018`).
- **Fig. S2D (equilibrium aspect ratio $L/r$ vs anchoring strength $\omega$)**: Digitized all 3 curves at 1200 dpi (`cal_res = (0.22, 0.16) px`, log-log axes):
  - `figS2d_numerical_min` (`d_figS2d_numerical_min`): 14 black solid circles (`omega = 0.1..4.0`, `L/r = 1.116..4.032`), reconstructing high-$\omega$ dot centers from the unobstructed top crown + measured isolated-dot radius where the green asymptotic curve overlaps the bottom rim.
  - `figS2d_scaling_small_w` (`d_figS2d_scaling_small_w`): 7 points along the red small-$\omega$ scaling curve `1 + omega` (`omega = 0.08..0.5`, `L/r = 1.080..1.500`).
  - `figS2d_scaling_large_w` (`d_figS2d_scaling_large_w`): 7 points along the green large-$\omega$ scaling curve `2*omega^(1/2)` (`omega = 1.2..4.0`, `L/r = 2.193..3.996`).
  - Inset colored droplet schematics left out.
  - **Model comparisons**: `cmp_figS2d_numerical_min` (`reproduced`, 14/14, max error `0.016`), `cmp_figS2d_scaling_small_w` (`reproduced`, 7/7, max error `0.0012`), `cmp_figS2d_scaling_large_w` (`reproduced`, 7/7, max error `0.0038`), `cmp_figS2d_large_omega_eq` (`reproduced`, 7/7, max error `0.0038`).
- **Fig. S3**: Left out — schematic geometry of cylindrical bundle capillary/pearling instability ($\rho(z)$, $R_0 - u^2/2$, $\lambda = 2\pi R_0 / k$, axial flow $v_z$); no quantitative axes.

## Pages 11–16 (printed pp. SI-5–SI-10) — Table S1 and Movies S1–S7
- **Table S1 (p. 11 / SI-5)**: Already transcribed in `data.yaml`; no figure panels.
- **Movie S1 (p. 11), Movies S2–S3 (p. 12), Movie S4 (p. 13), Movie S5 (p. 14), Movie S6 (p. 15), Movie S7 (p. 16)**: Left out — fluorescence microscopy video still frames and captions without quantitative plot axes.
