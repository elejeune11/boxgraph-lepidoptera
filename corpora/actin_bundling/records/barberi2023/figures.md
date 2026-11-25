# Figure Digitization Log: `barberi2023`

Source: Luca Barberi & Karsten Kruse, *Localized States in Active Fluids*, Phys. Rev. Lett. **131**, 238401 (2023) + 12-page Supplementary Material (PDF pp. 1–18).
Date: 2026-11-25.

All 52 calibrated figure readings (`806` total digitized data points across 10 quantitative figures / 20 panels) were extracted directly from high-resolution crops of `corpus/actin_bundling/2023_Barberi.pdf` via self-calibrating `mode: code` scripts in `figures/` (`_common.py`, `_panels.py`, and per-reading entrypoints). Every single overlay was generated under `figwork/overlays/` and visually verified on compact contact sheets (`figwork/verify_sheet_1.png`–`verify_sheet_6.png`, `figwork/verify_fixes.png`, `figwork/verify_aligned.png`); 100% of crosses sit on their intended curves or markers, and all readings are saved in `readings.json`.

## Summary of Key Quantitative Verdicts Settled by Digitization

1. **Main Text Fig. 1(b)–(g) vs. SM Table I (`q_sm_table1_W_Kobs_vs_fig1` settled):** All 10 curves of Fig. 1(a)–(g) (`262/262` points) are `reproduced` by converged 1D PDE solutions of Eqs. (SM1)–(SM5). Extracting FWHM width `W` directly from the calibrated Fig. 1(b)–(g) curves (`profile_fwhm_from_curve`) yields `W_fig = 7.36 (LP2), 6.79 (LP3), 11.00 (LP4), 10.82 (LP5), 10.36 (LP6), 11.41 (LP7)`, matching our PDE (`7.27, 6.78, 10.95, 10.83, 10.28, 11.38`) within `0.01–0.09` while contradicting SM Table I (`7.55, 7.92, 11.17, 11.60, 11.47, 12.58`; 6/6 `not reproduced`, off by up to `-1.17` / `-14%` on `LP3, LP5, LP6, LP7`). Together with Fig. 2, SM Fig. 4, and SM Fig. 5 (all matching our PDE), this proves SM Table I's widths and wavenumbers do not match the states plotted in Fig. 1.
2. **Main Text Fig. 3(c) spatial-dynamics eigenvalues (`q_fig3c_spatial_eigenvalues_mismatch` settled):** Calibrated broken-axis digitization measures the plotted markers at `Re = {0.001, 1.000, 30.527, 100.020}, Im = 0` (spike `×`) and `Re = {0.001, 1.000, 56.316, 56.316}, Im = {0, 0, -18.470, +18.406}` (`LP0` `+`), matching the decoupled 3x3 reaction-diffusion sub-block eigenvalues (`{30.5, 100}` and `{56.336 ± 18.464i}`) to `0.02–0.15%`, and returning `not reproduced` (`max|diff| = 92.41, 51.65, 14.56`) against the true 8x8 Jacobian of Eqs. (SM24)–(SM31) (`{0, 1.308, 2.862, 7.580}` and `{0, 0.670, 4.640 ± 3.891i}`).
3. **SM Fig. 1 stability boundary (`q_sm_fig1_long_wave_type_iis` refined):** Digitized `Omega_c(Z)` sits `+0.17` to `+0.92` above the continuum `K -> 0+` Type `II_s` boundary (`Omega_c = 8.39 [8.26, 8.53]` at `Z = 15` vs continuum `8.145`, `med|diff| = 0.375`), whereas imposing the finite periodic domain floor `K >= 2*pi/L = 0.2` yields `Omega_c = 8.32` at `Z = 15` and drops `med|diff|` nearly 3x to `0.134` (within the figure's `~0.2–0.3` scan-grid stair-stepping), consistent with Fig. 2's pink `H` branch extending to `Omega = 8.5`.
4. **SM Fig. 3 inset spike convergence (`q_audit_spike_grid_resolution` confirmed):** All 6 inset dots (`0.65272, 0.66710, 0.67374, 0.67694, 0.67850, 0.67928 ± 0.00005` at `L x 1..32`) sit uniformly `+0.0023..+0.0035` (`0.35–0.53%`) above cell-centered finite-volume `C(0.5*dx)` at fixed `Delta X = L/1024` (`0.64924..0.67697`, `not reproduced` at `±0.00005` precision) and `~0.012–0.014` below the continuum peak (`0.6644..~0.693`), confirming finite-difference grid under-resolution at `Delta X = L/1024`.

---

## Panel-by-Panel Log (All Figures Across PDF Pages 1–18)

### Main Text

- **Fig. 1(a) (PDF p. 2 / `238401-2`) — READ:**
  - Reading IDs: `fig1a_C` (solid black `C(X)`, 21 pts), `fig1a_Na` (dashed black `N_a(X)`, 21 pts), `fig1a_Ni` (dotted black `N_i(X)`, 21 pts), `fig1a_V` (pink `V(X)`, 21 pts).
  - Calibration: 1200 dpi crop `[5160, 740, 7300, 3860]`; 2 outward bottom ticks (`-L/2, +L/2`), 2 left ticks (`0, 3`), 2 pink right ticks (`-0.7, 0.7`); 2-tick exact linear calibration (0.00 px misfit). Curves separated by RGB color (pink vs dark) and connected-component span (`>=160 px` solid, `15–55 px` dashed with wing fallback where `N_a = C`, `3–14 px` dotted).
  - Comparisons: `cmp_fig1a_C`, `cmp_fig1a_Na`, `cmp_fig1a_Ni`, `cmp_fig1a_V` vs `pde_1d_profile_sm1_5` (`LP0`, `Omega = 6`): all 4 **`reproduced`** (`84/84` pts; `med|diff| = 0.0007, 0.0055, 0.022, 0.0015`).

- **Fig. 1(b)–(g) (PDF p. 2 / `238401-2`) — READ:**
  - Reading IDs: `fig1b_LP2_C` (27 pts), `fig1c_LP3_C` (27 pts), `fig1d_LP4_C` (31 pts), `fig1e_LP5_C` (31 pts), `fig1f_LP6_C` (31 pts), `fig1g_LP7_C` (31 pts).
  - Calibration: 1800 dpi crops per panel (`1106 px` x-span, `1101–1128 px` y-span); each panel individually calibrated from its own detected outward ticks (`-L/2, +L/2` and `0, 1`; note top tick `1` shifts panel-by-panel relative to the top frame due to per-panel y-padding). Curve isolated via largest connected component spanning the frame (eliminating `(b)`–`(g)` and `LP_i` text).
  - Comparisons: `cmp_fig1b_LP2_C` through `cmp_fig1g_LP7_C` vs `pde_1d_profile_sm1_5`: all 6 **`reproduced`** (`178/178` pts; `med|diff| = 0.0058–0.017`). Extracted FWHMs `cmp_W_fig1_vs_table1_LP2..LP7` are all 6 **`not reproduced`** against SM Table I while agreeing with our PDE `W` to `0.01–0.09`.

- **Fig. 2 (PDF p. 3 / `238401-3`) — READ:**
  - Reading IDs: `fig2_H` (pink dots, 13 pts), `fig2_LP0` (left green branch, 13 pts), `fig2_LP2` (blue dots, 4 pts), `fig2_LP3` (green dots, 6 pts), `fig2_LP4` (blue dots, 7 pts), `fig2_LP5` (green dots, 9 pts), `fig2_LP6` (blue dots, 12 pts), `fig2_LP7` (green dots, 8 pts) — 72 total points across all 8 plotted branches.
  - Calibration: 900 dpi crop `[615, 555, 3795, 2895]`; bottom ticks (`6, 19`, 2489 px span), left ticks (`0, 0.12`, 1998.5 px span). Circular markers (`radius ≈ 22.8 px`) separated from connecting lines by RGB color + 2x downsampled Euclidean distance transform with sub-pixel centroid refinement and slope-propagated horizontal uncertainty.
  - Comparisons: `cmp_fig2_H`, `cmp_fig2_LP2`, `cmp_fig2_LP3` **`reproduced`** (`23/23` pts); `cmp_fig2_LP0`, `cmp_fig2_LP4`, `cmp_fig2_LP5`, `cmp_fig2_LP6`, `cmp_fig2_LP7` **`close`** (`37/49` reproduced, `12/49` close, 0 not reproduced; `med|diff| = 0.00012–0.00029`, `max|diff| <= 0.00059` on a `0.12` y-scale).

- **Fig. 3(a) (PDF p. 3 / `238401-3`) — READ:**
  - Reading IDs: `fig3a_C` (solid black, 19 pts, avoiding pink-occluded `X = 0` column via `-0.06` and `+0.2`), `fig3a_Na` (dashed black, 14 pts on unobstructed `|X| >= 0.8`), `fig3a_Ni` (dotted black, 19 pts), `fig3a_V` (pink, 19 pts).
  - Calibration: 1200 dpi crop `[5200, 740, 7420, 2300]`; bottom ticks (`-L/2, +L/2`), left ticks (`0, 1.2`), pink right ticks (`-0.2, 0.2`).
  - Comparisons: `cmp_fig3a_C`, `cmp_fig3a_Na`, `cmp_fig3a_Ni`, `cmp_fig3a_V` vs `pde_1d_profile_sm1_5` (`spike`): all 4 **`reproduced`** (`71/71` pts; `med|diff| = 0.0002–0.0011`).

- **Fig. 3(b) (PDF p. 3 / `238401-3`) — READ:**
  - Reading IDs: `fig3b_spike_upper` (solid upper arc, 17 pts), `fig3b_spike_lower` (solid lower arc, 17 pts), `fig3b_LP0_upper` (dashed upper arc `C >= 0.58`, 11 pts), `fig3b_LP0_lower` (dashed lower arc `C <= 0.58`, 11 pts).
  - Calibration: 1200 dpi crop `[7420, 740, 9420, 2300]`; bottom ticks (`-4, 4`), left ticks (`0.2, 0.8`). Solid spike loop isolated by span `>= 250 px`, dashed `LP0` loop by `15–60 px` dashes split vertically at `C = 0.58`.
  - Comparisons: `cmp_fig3b_spike_upper`, `cmp_fig3b_spike_lower`, `cmp_fig3b_LP0_upper`, `cmp_fig3b_LP0_lower` vs `pde_1d_phase_orbit_sm1_5`: all 4 **`reproduced`** (`56/56` pts).

- **Fig. 3(c) (PDF p. 3 / `238401-3`) — READ:**
  - Reading IDs: `fig3c_spike_re` (4 pts), `fig3c_spike_im` (4 pts), `fig3c_LP0_re` (4 pts), `fig3c_LP0_im` (4 pts), indexed `1..4` ordered by `(Re, Im)`.
  - Calibration: 1200 dpi broken-axis crop `[5200, 2285, 9420, 3860]`; 5 left y-ticks (`-20, -10, 0, 10, 20`, misfit `0.25 px`), 2 left x-ticks (`0, 1`), 4 right x-ticks (`40, 60, 80, 100`, misfit `0.00 px`). `+` vs `×` markers distinguished automatically by horizontal/vertical arm ink density and overlapping `*` area.
  - Comparisons: vs `spatial_dynamics_spectrum_sm24` (exact 8x8 Jacobian of Eqs. SM24–SM31): `cmp_fig3c_spike_im` **`reproduced`** (`4/4` real), while `cmp_fig3c_spike_re` (`max|diff| = 92.41`), `cmp_fig3c_LP0_re` (`max|diff| = 51.65`), and `cmp_fig3c_LP0_im` (`max|diff| = 14.56`) are **`not reproduced`**.

- **Fig. 4(a)–(d) (PDF p. 4 / `238401-4`) — LEFT OUT:**
  - Reason: 2D spatial colormap snapshots of `C(X, Y)` (colorbar `0.2` to `1.1`) overlaid with 2D velocity arrowheads in `(a), (c)`; no 1D curve, scatter points, or bar values to digitize. (2D PDE verification is already covered by `run_pde_2d_fig4a` and `run_pde_2d_fig4a_printed`.)

### Supplementary Material

- **SM Fig. 1 (PDF p. 11 / `SM-5`) — READ:**
  - Reading IDs: `smfig1_boundary_Io` (dashed Type `I_o` boundary `Z in [0.6, 2.6]`, 6 pts), `smfig1_boundary_Is` (solid boundary `Z in [3.5, 24.2]`, 17 pts avoiding dot-on-curve columns), `smfig1_LS_onset` (lower existence envelope of the 80 circular black `LS` dots across `Z = 3, 6, ..., 30`, 10 pts, `origin: run`).
  - Calibration: 600 dpi crop `[560, 540, 4540, 3560]`; bottom ticks (`0, 30`, 3404 px span), left ticks (`0, 30`, 2501.5 px span). Circular `LS` dots (`radius ≈ 20 px`) detected via 2x downsampled EDT and masked before reading the stability curve.
  - Comparisons: `cmp_smfig1_boundary_Io` **`close`** (`3/6` rep, `3/6` close; `med|diff| = 0.13`, `max|diff| = 0.19`); `cmp_smfig1_boundary_Is_cont` (`k_min = 0.02`, `med|diff| = 0.375`) and `cmp_smfig1_boundary_Is_finiteL` (`k_min = 0.20 = 2*pi/L`, `med|diff| = 0.134`) both **`reproduced`** after input rounding propagation, with finite-domain `K_1 = 0.2` matching 3x closer to the plotted curve.

- **SM Fig. 2 (PDF p. 12 / `SM-6`) — LEFT OUT:**
  - Reason: 1D space-time kymograph `C(X, T)` (`T in [0, 2000]`, colorbar `0.3` to `0.7`) of a passive-fluid travelling wave; kymographs have no 1D curve or discrete point series to read (excluded per guide; PDE travelling-wave verification is in `run_pde_wave_smfig2`).

- **SM Fig. 3 Main Panel (PDF p. 15 / `SM-9`) — READ:**
  - Reading IDs: `smfig3_main_C_L32` (outermost gray `L x 32` spike profile over `X in [-L/8, L/8]`, 19 pts; note individual `L x 2..32` curves overlap completely outside the immediate peak tip), `smfig3_main_H` (pink dashed homogeneous steady state `H`, 7 pts).
  - Calibration: 600 dpi crop `[560, 560, 4500, 3540]`; bottom ticks (`-L/8, +L/8`, 3577 px span), left ticks (`0.2, 0.7`, 2273 px span). Direct pixel check shows the spike peak sits at `X = -0.0165 ≈ -Delta X / 2`.
  - Comparisons: `cmp_smfig3_main_C_L32` (`19/19`, `med|diff| = 0.0014`) and `cmp_smfig3_main_H` (`7/7`, `med|diff| = 0.00005`) both **`reproduced`**.

- **SM Fig. 3 Inset (PDF p. 15 / `SM-9`) — READ:**
  - Reading ID: `smfig3_inset_Cmax` (all 6 black circular dots `Max. C` vs `System size (x L)` at `1, 2, 4, 8, 16, 32`).
  - Calibration: 1200 dpi crop `[5720, 1400, 8740, 3720]`; 16 evenly spaced bottom ticks (`2, 4, ..., 32`, misfit `0.00 px`) and 4 left ticks (`0.65, 0.66, 0.67, 0.68`, misfit `0.00 px`). Read values: `0.65272, 0.66710, 0.67374, 0.67694, 0.67850, 0.67928` (`±0.00005`).
  - Comparison: `cmp_smfig3_inset_Cmax` vs `pde_1d_spike_scaling_smfig3` (`Delta X = L/1024`): **`not reproduced`** (`0/6` inside `±0.00005`, uniform cell-center vs vertex-center offset `-0.0023` to `-0.0035`, `0.35–0.53%`).

- **SM Fig. 4 Left (`LP0`) & Right (`LP3`) (PDF p. 16 / `SM-10`) — READ:**
  - Reading IDs: `smfig4_LP0_neg_dX_Jc` (solid black, 17 pts on unobstructed `|X| <= 5`), `smfig4_LP0_Rc` (dotted black, 17 pts on `|X| <= 5`), `smfig4_LP0_V` (pink, 21 pts on `[-12, 12]`); `smfig4_LP3_neg_dX_Jc` (solid black, 19 pts on `|X| <= 5.5`), `smfig4_LP3_Rc` (dotted black, 19 pts on `|X| <= 5.5`), `smfig4_LP3_V` (pink, 23 pts on `[-12, 12]`). (Outside `|X| > 6`, pink `V = 0` is painted over `-d_X J_c = R_c = 0`.)
  - Calibration: 600 dpi crops (`[590, 580, 2530, 3540]` left, `[2530, 580, 4490, 3540]` right); bottom ticks (`-L/2, +L/2`), left ticks (`-1, 1`), pink right ticks (`-0.7, 0.7`).
  - Comparisons: all 6 (`cmp_smfig4_LP0_neg_dX_Jc`, `Rc`, `V`, `cmp_smfig4_LP3_neg_dX_Jc`, `Rc`, `V`) **`reproduced`** (`116/116` pts).

- **SM Fig. 5(a)–(d) (PDF p. 17 / `SM-11`) — READ:**
  - Reading IDs: Panel (a) `smfig5a_C` (15 pts), `smfig5a_Na` (15 pts), `smfig5a_Ni` (13 pts on `[14.75, 16.25]` outside top-right legend), `smfig5a_V` (15 pts); Panel (b) `smfig5b_Jca` (pink, 15 pts), `smfig5b_Jcd` (gray, 15 pts); Panel (c) `smfig5c_Jaa` (pink, 15 pts), `smfig5c_Jad` (gray, 15 pts); Panel (d) `smfig5d_Jia` (pink, 15 pts), `smfig5d_Jid` (gray, 15 pts).
  - Calibration: 600 dpi crop `[640, 550, 4520, 2490]` for (a) (`X in [14.7, 16.6]`, left `[0.2, 2]`, right `[-0.2, 0.2]`); 1200 dpi crops for (b) `[-1, 0.5]`, (c) `[-2, 1]`, (d) `[-2, 5]`. Legend swatches and panel letters masked surgically while preserving faded curves visible behind the translucent legend boxes.
  - Comparisons: all 10 (`cmp_smfig5a_C`, `Na`, `Ni`, `V`, `cmp_smfig5b_Jca`, `Jcd`, `cmp_smfig5c_Jaa`, `Jad`, `cmp_smfig5d_Jia`, `Jid`) **`reproduced`** (`148/148` pts; `med|diff| = 0.00038` on `C`, `0.00081` on `V`, `0.0010` on `-d_X J_c^d`) by `LP3` at `Omega = 10` aligned to the symmetric central peak (`X_0 = 15.661` from `ref_C`).

- **SM Fig. 6 Left & Right (PDF p. 18 / `SM-12`) — LEFT OUT:**
  - Reason: 2D zoomed colormap snapshots (`C in [0.2, 1.1]`) with 2D black velocity arrowheads corresponding to Fig. 4(b) and 4(d); no 1D curve or scatter values to digitize.
