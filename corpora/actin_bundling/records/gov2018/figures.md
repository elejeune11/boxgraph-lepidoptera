# Figure Readings — `gov2018`

Complete inventory and calibrated `mode: code` digitization of every figure panel across all 13 pages (Figs. 1–6) of N. S. Gov, *Phil. Trans. R. Soc. B* **373**: 20170115 (2018) (`corpus/actin_bundling/2018_Gov.pdf`). All quantitative 2D panels were rendered at **1200 dpi** (**1800 dpi** for the compact Fig. 4(d) inset), axis frames and inward ticks were detected directly from the rendered pixels via `figures/_calib.py`, and every check overlay PNG in `figwork/` (`boxgraph figure check`) was visually inspected to confirm all magenta crosses lie squarely on their target curves.

---

## Pages Without Figures

- **Page 1:** Title, abstract, and Section 1 (Introduction); no figures.
- **Page 3:** Section 2 equations (2.4)–(2.9) and prose; no figures.
- **Page 5:** Sections 3–4 prose; no figures.
- **Page 7:** Sections 4–5 prose and Eq. (5.1); no figures.
- **Pages 10–13:** Section 6 (Discussion) and References 1–132; no figures.

---

## Complete Figure Inventory (Figs. 1–6)

### Figure 1 (p. 2) — Left out
- **Fig. 1(a) — Left out:** Uncalibrated schematic of a localized membrane height deformation $h(x)$ vs. $x$ (arrows only; no ticks or numerical axis scales).
- **Fig. 1(b) — Left out:** Uncalibrated schematic of corresponding mean curvature $H(x)$ vs. $x$ (no ticks or numerical axis scales).
- **Fig. 1(c) — Left out:** Schematic cartoon illustrating lateral flow ($\mathbf{J}_{\text{curv}}$) and cytoplasmic adsorption of convex (blue) and concave (red) CMC; no axes.

### Figure 2 (p. 4) — Left out
- **Fig. 2(a) — Left out:** Uncalibrated schematic of positive feedback for convex CMC recruiting protrusive forces ($\lambda = 2\pi/q_{\max}$, symbolic time arrow; no numerical axes).
- **Fig. 2(b) — Left out:** Qualitative sketch of linear stability growth rate $\omega$ vs. wavevector $q$ showing only symbolic $q_{\max}$; axes have no tick marks and no numerical scale values.
- **Fig. 2(c) — Left out:** Uncalibrated 2D closed-contour simulation snapshots (dotted, dashed, solid shapes from Kabaso *et al.* [35]); no coordinate axes, ticks, or scale bar.
- **Fig. 2(d) — Left out:** Uncalibrated steady-state fan-shaped closed membrane contour from [35]; no axes, ticks, or scale bar.
- **Fig. 2(e) — Left out:** Fluorescence micrograph of actin in a fan-shaped cell (10 $\mu\text{m}$ scale bar only) with two overlaid white theoretical contours from [35]; no quantitative graph axes.
- **Fig. 2(f) — Left out:** HIV viral budding cartoons plus fluorescence micrograph and AFM height image (4.0 $\mu\text{m}$ scale bar only) from Gladnikoff *et al.* [58]; no quantitative graph axes.

### Figure 3 (p. 6) — Left out
- **Fig. 3(a) — Left out:** Uncalibrated schematic of circumferential arc-like CMC rings on a cylindrical tube ($z$, $h(z)$, $n(z)$); no numerical scales.
- **Fig. 3(b) — Left out:** Uncalibrated schematic of axially aligned CMC on a constricted tube of length $L$; no numerical scales.
- **Fig. 3(c)(i)–(iv) — Left out:** Four 3D oblique perspective MATLAB surface renderings of cylindrical ring coarsening at $t = 0.002, 0.1, 1, 2$ from Shlomovitz & Gov [88]. Although the 3D bounding boxes carry bare numbers (`0, 10, 20` axially and `-0.5, 0, 0.5` transversely), the axes have no labels or units, the perspective projection lacks orthogonal 2D axes, and the cyan/magenta surface shading has no colorbar; the four printed snapshot times are already recorded in `data.yaml` (`fig3c_time_i`–`iv`).
- **Fig. 3(d) — Left out:** Uncalibrated schematic of a cortical actin–myosin network squeezing a membrane tube along axis $Z$; no numerical scales.

### Figure 4 (p. 8) — Panels (a)–(b) left out; Panels (c)–(d) and Inset top read

- **Fig. 4(a) — Left out:** Uncalibrated schematic of convex CMC recruiting actin ($A$) and adsorbing myosin ($k_{\text{on}}, k_{\text{off}}, A^*$) from [116]; no numerical axes.
- **Fig. 4(b) — Left out:** Uncalibrated schematic of strong myosin retraction driving CMC aggregation onto dip shoulders and sideways wave propagation; no numerical axes.

- **Fig. 4(c) (p. 8) — Read (10 readings: `fig4c_omega_real_cyan`, `fig4c_omega_real_nonlocal`, `fig4c_omega_real_green`, `fig4c_omega_real_red`, `fig4c_omega_real_blue`, `fig4c_omega_imag_purple`, `fig4c_omega_imag_cyan`, `fig4c_omega_imag_nonlocal`, `fig4c_omega_imag_red`, `fig4c_omega_imag_green`):**
  - *Render:* `dpi: 1200`, `crop: [3360, 860, 6240, 3380]` (plot boxes 2184 × 876 px top, 2184 × 888 px bottom).
  - *Top panel $\omega'(q)$ (Hz) vs. $q$ ($\mu\text{m}^{-1}$):*
    - Calibration: 11 x-ticks (`0, 2, ..., 20`, misfit `3.17 px` = `0.029 um^-1`); 8 left-border y-ticks (`0.02, 0, -0.02, ..., -0.12 Hz`, misfit `4.27 px` = `0.00079 Hz`).
    - Auxiliary calibrated reference lines: horizontal dash-dotted detachment-rate line reads $\omega' = -0.09996\ \text{Hz}$ (confirming $k_{\text{off}} = 0.100\ \text{Hz}$); vertical dashed wave-onset line reads $q_\omega = 18.64\ \mu\text{m}^{-1}$ (top) / $18.68\ \mu\text{m}^{-1}$ (bottom).
    - `fig4c_omega_real_cyan` (`omega_real_cyan_fig4c`, 12 pts, $q \in [0.4, 19.5]$): solid light-blue curve (largest $A^*$); small-$q$ dip $-0.0076\ \text{Hz}$ at $q = 0.4$, broad wave-instability peak $+0.0264\ \text{Hz}$ at $q = 6$, zero crossing near $q = 14.8$, decreasing to $-0.0632\ \text{Hz}$ at $q = 19.5$.
    - `fig4c_omega_real_nonlocal` (`omega_real_nonlocal_fig4c`, 11 pts, $q \in [2, 19.5]$): dashed light-blue curve (non-local myosin kernel [116]), isolated below the solid cyan trace and sampled by Gaussian-weighted local quadratic regression; peaks near $+0.0144\ \text{Hz}$ at $q = 2$ and merges onto the red short-wavelength branch for $q \ge 14$.
    - `fig4c_omega_real_green` (`omega_real_green_fig4c`, 11 pts, $q \in [3, 19.6]$): solid green curve (smallest nonzero $A^*$); stationary-protrusion peak $+0.0108\ \text{Hz}$ at $q = 9$, drops to $-0.0908\ \text{Hz}$ at $q = 18.2$, then flattens near $-k_{\text{off}} = -0.101\ \text{Hz}$ for $q > q_\omega$.
    - `fig4c_omega_real_red` (`omega_real_red_fig4c`, 12 pts, $q \in [2, 19.5]$): solid red curve (intermediate $A^*$); decays smoothly to the damped-wave bifurcation kink at $(q \approx 9.34, \omega' = -0.0449\ \text{Hz})$ and down to $-0.1123\ \text{Hz}$ at $q = 19.5$.
    - `fig4c_omega_real_blue` (`omega_real_blue_fig4c`, 4 pts, $q \in [18.8, 19.4]$): solid dark-blue curve ($A^* = 0$) read on its visible branch below $-k_{\text{off}}$ ($-0.1049$ to $-0.1254\ \text{Hz}$); for $q \le 18.6$ it is obscured directly behind the green curve.
  - *Bottom panel $\omega''(q)$ (Hz) vs. $q$ ($\mu\text{m}^{-1}$):*
    - Calibration: 11 x-ticks (`0, 2, ..., 20`, misfit `1.47 px` = `0.013 um^-1`); left y-axis (`0.15, 0.10, 0.05, 0.00 Hz`, exact 282.0 px step, misfit `0.00 px`, zero at row `2095.0`); right y-axis (`5, 4, 3, 2, 1, 0 Hz`, exact 168.0 px step, misfit `0.00 px`, zero at row `2089.0`).
    - *Right-axis resolution (`q_fig4c_right_axis_assignment`):* Gov (2018)'s caption omits the label and curve mapping of the right-hand `0–5` axis. Subpixel zero-intercept alignment (row `2089` vs. `2095`), mutual tangency at $q \to 0$, and 3×3 cubic dispersion consistency with $k_{\text{off}} = 0.100\ \text{Hz}$ and $\omega'_{\max} = +0.026\ \text{Hz}$ show that purple, solid light blue, and dashed light blue belong to the right `0–5 Hz` scale, whereas red and green sit on the left `0–0.15 Hz` scale (dark blue $A^* = 0$ is identically zero throughout).
    - `fig4c_omega_imag_purple` (`omega_imag_purple_fig4c`, 9 pts, $q \in [4, 19.5]$, right axis): analytical linear approximation from $0.999\ \text{Hz}$ ($q = 4$) to $4.853\ \text{Hz}$ ($q = 19.5$), slope $c = 0.2499\ \mu\text{m}\cdot\text{Hz}$.
    - `fig4c_omega_imag_cyan` (`omega_imag_cyan_fig4c`, 10 pts, $q \in [2, 19.5]$, right axis): largest-$A^*$ solid light-blue curve from $0.428\ \text{Hz}$ ($q = 2$) to $4.356\ \text{Hz}$ ($q = 19.5$), bending softly below the purple tangent at high $q$.
    - `fig4c_omega_imag_nonlocal` (`omega_imag_nonlocal_fig4c`, 9 pts, $q \in [5, 19]$, right axis): non-local kernel dashed light-blue curve peaking at $1.057\ \text{Hz}$ near $q = 7.5$ and decaying to $0.088\ \text{Hz}$ at $q = 19$.
    - `fig4c_omega_imag_red` (`omega_imag_red_fig4c`, 11 pts, $q \in [2, 19.5]$, left axis): $0.000\ \text{Hz}$ for $q \le 9$, bifurcating sharply above $q \approx 9.34$ to $0.1064\ \text{Hz}$ at $q = 19.5$.
    - `fig4c_omega_imag_green` (`omega_imag_green_fig4c`, 8 pts, $q \in [12, 19.5]$, left axis): $0.000\ \text{Hz}$ outside the narrow wave-onset arch straddling $q_\omega = 18.66\ \mu\text{m}^{-1}$, peaking at $0.0054\ \text{Hz}$ ($q = 18.5, 18.8$).
  - *Model comparisons:*
    - `cmp_fig4c_purple_acoustic` (`purpose: reproduction`, swept on `myosin_acoustic_dispersion`, `prefactor: unknown`): **reproduced** (9/9 points reproduced; implied prefactor `1.443` relative to baseline `myosin_c_ana_wave = 0.1732 um/s`, giving exact linear slope $c = 0.2499\ \mu\text{m}\cdot\text{Hz}$).
    - `cmp_fig4c_cyan_acoustic` (`purpose: approximation`, swept on `myosin_acoustic_dispersion`, `prefactor: unknown`): **not reproduced** (4 reproduced, 4 close, 2 not reproduced at $q = 10, 19.5$; implied prefactor `1.346`, $c_{\text{fit}} = 0.233\ \mu\text{m}\cdot\text{Hz}$), quantifying the $-10.6\%$ downward bending of $\omega''(q)$ below the small-$q$ linear tangent at $q = 19.5\ \mu\text{m}^{-1}$ while confirming Gov's p. 8 description that it is "almost linear" at longer wavelengths.
    - The remaining Fig. 4(c) curves cannot be compared against `myosin_wave_dispersion` without unstated simulation parameters (`q_unstated_simulation_parameters`) or the omitted non-local myosin kernel (`q_nonlocal_myosin_kernel_unstated`).

- **Fig. 4(d) main (p. 8) — Read (10 readings: `fig4d_h_t0`–`t3`, `fig4d_m_t1`–`t3`, `fig4d_n_t1`–`t3`):**
  - *Render:* `dpi: 1200`, `crop: [6300, 820, 8980, 3520]` (box width 1953 px; heights 1266, 528, 534 px). The inset box inside top panel $h(x)$ is auto-detected and masked out.
  - *Calibration:* Shared x-axis (`0, 0.5, ..., 2.5 um`, misfit `1.43 px` = `0.0018 um`); $h(x)$ left y-axis (`0` to `-0.6 um`, misfit `2.25 px` = `0.0012 um`); $h(x)$ right y-axis for $t = 0$ (`0` anchored to zero line, `-0.5` and `-1` from right-side label centroids, misfit `4.50 px` = `0.0039 um`); $m(x)$ y-axis (`0.5, 0, -0.5, -1.0 um^-3`, misfit `0.90 px` = `0.0025 um^-3`); $n(x)$ y-axis (`1.0, 0.5, 0, -0.5, -1.0, -1.5 um^-3`, misfit `0.95 px` = `0.0045 um^-3`).
  - *Curves:*
    - `fig4d_h_t0` (`h_t0_fig4d`, 8 pts, $x \in [0.05, 0.8]$): initial Gaussian height dip from $-0.987\ \mu\text{m}$ ($x = 0.05$) to $-0.003\ \mu\text{m}$ ($x = 0.8$) on right-hand scale.
    - `fig4d_h_t1`, `fig4d_h_t2`, `fig4d_h_t3` (`h_t1_fig4d`–`h_t3_fig4d`, 10 pts each): $h(0.05) = -0.617, -0.494, -0.396\ \mu\text{m}$ with propagating leading crests $+0.040, +0.040, +0.037\ \mu\text{m}$ near $x = 1.1, 1.35, 1.6\ \mu\text{m}$.
    - `fig4d_m_t1`, `fig4d_m_t2`, `fig4d_m_t3` (`m_t1_fig4d`–`m_t3_fig4d`, 7/8/9 pts): myosin density peaks $+0.347, +0.433, +0.433\ \mu\text{m}^{-3}$ near $x = 0.58, 0.68, 0.76\ \mu\text{m}$. ($t = 0$ is flat $m \equiv 0$, left out as trivial initial baseline.)
    - `fig4d_n_t1`, `fig4d_n_t2`, `fig4d_n_t3` (`n_t1_fig4d`–`n_t3_fig4d`, 9/8/10 pts): actin/CMC density peaks $+0.946, +0.552, +0.397\ \mu\text{m}^{-3}$ near $x = 0.78, 1.02, 1.22\ \mu\text{m}$ — confirming quantitative actin lead over myosin of $+0.20, +0.34, +0.46\ \mu\text{m}$ at $t = 1, 2, 3\ \text{s}$. ($t = 0$ is flat $n \equiv 0$, left out.)
  - *Model comparisons:* None directly, because simulation parameters and the dimensionless-to-$\mu\text{m}^{-3}$ density conversion factor are unstated in Gov (2018) (`q_unstated_simulation_parameters`, `q_n_unit_shift_fig4d`).

- **Fig. 4(d) inset (p. 8) — Top $h(x)$ read (`fig4d_inset_h_t60`, `fig4d_inset_h_t90`, `fig4d_inset_h_t120`); middle $m(x)$ and bottom $n(x)$ left out:**
  - *Render:* `dpi: 1800`, `crop: [10830, 1500, 13050, 3150]` (inset $h(x)$ plot box 1651 × 706 px).
  - *Calibration:* 6 top x-ticks (`0, 5, ..., 25 um`, misfit `2.10 px` = `0.032 um`); 4 right y-ticks (`0.05, 0, -0.05, -0.10 um`, misfit `1.80 px` = `0.00063 um`).
  - *Curves (9 pts each):*
    - `fig4d_inset_h_t60` (`h_inset_t60_fig4d`): trough $-0.130\ \mu\text{m}$ at $x = 7.4\ \mu\text{m}$, crest $+0.045\ \mu\text{m}$ at $x = 9.9\ \mu\text{m}$.
    - `fig4d_inset_h_t90` (`h_inset_t90_fig4d`): trough $-0.105\ \mu\text{m}$ at $x = 11.0\ \mu\text{m}$, crest $+0.046\ \mu\text{m}$ at $x = 13.7\ \mu\text{m}$.
    - `fig4d_inset_h_t120` (`h_inset_t120_fig4d`): trough $-0.088\ \mu\text{m}$ at $x = 14.8\ \mu\text{m}$, crest $+0.053\ \mu\text{m}$ at $x = 17.8\ \mu\text{m}$ (yielding long-time wave propagation speed $0.132\ \mu\text{m/s}$ for crests and $0.123\ \mu\text{m/s}$ for troughs).
  - *Inset middle $m(x)$ and bottom $n(x)$ — Left out:* Neither inset subpanel carries numeric y-axis tick labels or scale values (uncalibrated vertical axes).

### Figure 5 (p. 9) — Panels (a)–(c) left out; Panel (d) read

- **Fig. 5(a) — Left out:** Uncalibrated schematic of diffusing convex ($-$) and concave ($+$) CMC with active forces $A_-, A_+$; no numerical axes.
- **Fig. 5(b) — Left out:** Uncalibrated schematic of cytoplasmic adsorption/desorption ($k_{\text{on}}, k_{\text{off}}$) of convex/concave CMC; no numerical axes.
- **Fig. 5(c) — Left out:** Uncalibrated feedback-loop diagram linking $h$, $n_-$, and $n_+$; no axes.

- **Fig. 5(d) (p. 9) — Read (3 readings: `fig5d_h`, `fig5d_n_convex`, `fig5d_n_concave`):**
  - *Render:* `dpi: 1200`, `crop: [5480, 800, 9120, 2960]` (three stacked plot boxes 2960 × 581 px each). Sampled at 13 points across $x \in [-9.5, 4.7]\ \mu\text{m}$ covering nearly three full wavelengths while avoiding the right-side text/icon occlusions on $x \in [5.2, 9.2]\ \mu\text{m}$.
  - *Calibration:* 11 shared x-ticks (`-10, -8, ..., +10 um`, misfit `2.23 px` = `0.015 um`); $h(x)$ y-axis (`1.0, 0.5, 0, -0.5, -1.0 um`, misfit `2.20 px` = `0.0076 um`); $n_-(x)$ y-axis (`0.815, 0.810, 0.805, 0.800, 0.795, 0.790`, misfit `2.17 px` = `0.000093`); $n_+(x)$ y-axis (`0.16, 0.15, 0.14, 0.13, 0.12, 0.11`, misfit `1.96 px` = `0.00017`).
  - *Curves:*
    - `fig5d_h` (`h_wave_fig5d`, 13 pts): steady-state leftward wave ($\lambda = 5.0\ \mu\text{m}$, $V = 0.67\ \mu\text{m/min}$ printed), oscillating between trough $-0.497\ \mu\text{m}$ and crest $+0.863\ \mu\text{m}$.
    - `fig5d_n_convex` (`n_convex_wave_fig5d`, 13 pts): diffusing convex CMC $n_-(x)$ oscillating between $0.7924$ and $0.8079$, leading $h(x)$ crests slightly by $\sim 0.15\ \mu\text{m}$ ($\sim 11^\circ$) in the leftward propagation direction.
    - `fig5d_n_concave` (`n_concave_wave_fig5d`, 13 pts): adsorbing concave CMC $n_+(x)$ oscillating between $0.119$ and $0.152$, trailing $h(x)$ troughs by $\sim 0.7\ \mu\text{m}$ ($\sim 50^\circ$ off exact $180^\circ$ anti-phase; logged in `q_fig5d_phase_lag_vs_caption`).
  - *Model comparisons:* None possible because the two-curvature PDEs and parameters from Peleg *et al.* [41] are omitted from Gov (2018) (`two_curvature_wave_model_peleg2011` is `runnable: false`, `q_two_curvature_model_unstated`).

### Figure 6 (p. 9) — Left out
- **Fig. 6 (main and subcases i–iii) — Left out:** Uncalibrated mechano-chemical feedback network diagram linking calcium concentration $C$, membrane height $h$, and CMC density $n$ from Veksler & Gov [123]; no axes or quantitative values.
