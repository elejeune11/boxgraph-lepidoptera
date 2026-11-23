# Figure Readings — `bois2011`

Calibrated `mode: code` digitization of all quantitative 2D panels in Bois, Jülicher & Grill, *Phys. Rev. Lett.* **106**, 028103 (2011) (`corpus/actin_bundling/2011_Bois.pdf`, pages 1–4). All panels were rendered at **1200 dpi** (plot boxes 1054–1617 px per side), axis frames and inward ticks were detected directly from the rendered pixels via `figures/_calib.py`, and every overlay PNG in `figwork/` was visually inspected after `boxgraph figure check` to confirm all magenta crosses lie on their intended curves.

- **Page 1 (no figures):** Text and Eqs. (1)–(5) only; no figure panels on page 1.

- **Fig. 1(a) (p. 2) — Read (`fig1a_c`, `fig1a_v`, `fig1a_sigma`, `fig1a_zdm`):**
  - *Curves:* 1-peak nonhomogeneous steady state ($\mathrm{Pe} = 25$, $L/\ell = 2\pi$, $c_0 = 1$, periodic BCs): top subpanel concentration $c(x/\ell)$ (solid, left axis, 11 points) and flow velocity $v/U(x/\ell)$ (dashed, right axis, 11 points); bottom subpanel total stress $\sigma/(\zeta\Delta\mu)_0$ (solid, 10 points) and active stress $\zeta\Delta\mu/(\zeta\Delta\mu)_0$ (dashed, 10 points).
  - *Calibration (1200 dpi):* Frame border + 6 interior x-ticks ($x/\ell = 0..6$, residual `0.00 px`); top left 4 y-ticks ($c = 2, 4, 6, 8$, residual `0.15 px`); top right 3 y-ticks ($v/U = -0.2, 0, 0.2$, residual `0.00 px`); bottom 5 y-ticks ($\text{stress} = 0..0.8$, residual `0.30 px`).
  - *Reproduction verdicts:* Evaluated via `single_species_steady_state_profile_eq9` (wrapper over existing `single_species_steady_state_eq9`, peak centered at $x/\ell = \pi$): `cmp_fig1a_c` **reproduced**, `cmp_fig1a_v` **reproduced**, `cmp_fig1a_sigma` **reproduced**, `cmp_fig1a_zdm` **reproduced** (4/4 reproduced; peak values $c_{\max} = 8.245$ vs. $8.271$, $(v/U)_{\max} = 0.250$ vs. $0.250$, $\hat{\sigma}_{\max} = 0.4262$ vs. $0.4262$, $(\zeta\Delta\mu)_{\max} = 0.8935$ vs. $0.8921$).

- **Fig. 1(b) (p. 2) — Read (`fig1b_c`, `fig1b_v`, `fig1b_sigma`, `fig1b_zdm`):**
  - *Curves:* 2-peak nonhomogeneous steady state ($\mathrm{Pe} = 25$, $L/\ell = 2\pi$, $c_0 = 1$, periodic BCs): top subpanel $c(x/\ell)$ (solid, 12 points) and $v/U(x/\ell)$ (dashed, 11 points); bottom subpanel $\sigma/(\zeta\Delta\mu)_0$ (solid, 11 points) and $\zeta\Delta\mu/(\zeta\Delta\mu)_0$ (dashed, 11 points).
  - *Calibration (1200 dpi):* Same tick structure as Fig. 1(a); residuals `(0.00, 0.15) px` for $c$, `(0.00, 0.00) px` for $v/U$, `(0.00, 0.30) px` for stresses.
  - *Reproduction verdicts:* Evaluated via `single_species_steady_state_profile_eq9` (peaks at $x/\ell = \pi/2, 3\pi/2$): `cmp_fig1b_c` **reproduced**, `cmp_fig1b_v` **reproduced**, `cmp_fig1b_sigma` **reproduced**, `cmp_fig1b_zdm` **reproduced** (4/4 reproduced; e.g. $c_{\max} = 3.093$ vs. $3.095$, $c_{\min} = 0.100$ vs. $0.099$, $\hat{\sigma}_{\max} = 0.4511$ vs. $0.4512$, $\hat{\sigma}_{\min} = 0.3135$ vs. $0.3137$).

- **Fig. 2(a) (p. 2) — Read (`fig2a_Pe10`, `fig2a_Pe7`, `fig2a_Pe3`):**
  - *Curves:* Dimensionless Hamiltonian potentials $\ell^2 V(u)$ vs. $u = \ln c$ for $c_0 = 1$: $\mathrm{Pe} = 10$ (solid, 13 points), $\mathrm{Pe} = 7$ (dashed, 12 points), and $\mathrm{Pe} = 3$ (dotted, 8 points). Legend and inset frames are automatically detected and masked; dashed/dotted curves use local Gaussian-weighted quadratic pixel regression to avoid slanted dash-tip bias. (The inset magnifies the same three curves on $u \in [-2, 2]$ around $(0,0)$ and adds no distinct curve.)
  - *Calibration (1200 dpi):* Left/right borders + 3 interior x-ticks ($u = -8, -4, 0, 4, 8$, residual `0.40 px`); 5 right-border y-ticks ($\ell^2 V = -10..10$, residual `0.20 px`).
  - *Reproduction verdicts:* Swept directly on existing model `single_species_potential_eq10` (using $\hat{\sigma}_0 = 1/2$ for homogeneous $\mathrm{Pe} = 3, 7$ and `sigma0_hat_Pe10_pred` from `run_ss_fig2a_Pe10` for $\mathrm{Pe} = 10$): `cmp_fig2a_Pe10` **reproduced**, `cmp_fig2a_Pe7` **reproduced**, `cmp_fig2a_Pe3` **reproduced** (3/3 reproduced).

- **Fig. 2(b) (p. 2) — Read (`fig2b_2peak`, `fig2b_1peak`):**
  - *Curves:* Dimensionless potentials $\ell^2 V(u)$ vs. $u = \ln c$ for the $\mathrm{Pe} = 25$ steady states of Fig. 1: 2-peak steady state (solid, 12 points) and 1-peak steady state (dashed, 10 points).
  - *Calibration (1200 dpi):* Borders + 2 interior x-ticks ($u = -10, -5, 0, 5$, residual `0.15 px`); 5 right-border y-ticks ($\ell^2 V = 0..40$, residual `0.20 px`).
  - *Reproduction verdicts:* Swept directly on existing model `single_species_potential_eq10` using `sigma0_hat_fig1b_pred` and `sigma0_hat_fig1a_pred`: `cmp_fig2b_2peak` **reproduced**, `cmp_fig2b_1peak` **reproduced** (2/2 reproduced).

- **Fig. 3(a) (p. 3) — Read (`fig3a_Pe25`, `fig3a_Pe10`, `fig3a_Pe7`, `fig3a_Pe3`):**
  - *Curves:* Single-species dispersion relations $\tau_D \lambda(k)$ vs. $k\ell$ for $c_0 = 1$, $f(c) = c/(1+c)$: $\mathrm{Pe} = 25$ (top, 13 points), $\mathrm{Pe} = 10$ (second, 9 points), $\mathrm{Pe} = 7$ (third, 8 points), and $\mathrm{Pe} = 3$ (bottom, 7 points). Vertical gray lines ($k\ell = 1, 2$) are excluded by intensity threshold (`gray < 95`), and the thin zero-line, arrow, and "increasing Pe" text are removed by morphological opening before continuous curve tracing.
  - *Calibration (1200 dpi):* Borders + 2 interior x-ticks ($k\ell = 0, 1, 2, 3$, residual `0.00 px`); 5 right-border y-ticks ($\tau_D\lambda = -2..2$, residual `0.20 px`).
  - *Reproduction verdicts:* Swept directly on existing model `single_species_dispersion_hill`: `cmp_fig3a_Pe25` **reproduced**, `cmp_fig3a_Pe10` **reproduced**, `cmp_fig3a_Pe7` **reproduced**, `cmp_fig3a_Pe3` **reproduced** (4/4 reproduced; e.g. $\tau_D\lambda(1) = 2.124$ vs. exact $2.125$ at $\mathrm{Pe} = 25$).

- **Fig. 3(b) (p. 3) — Read (`fig3b_stability_boundary`):**
  - *Curves:* Log-log linear stability boundary $\mathrm{Pe}/(1 + (\beta\pi\ell/L)^2)$ read horizontally at 11 log-spaced average concentrations $c_0 \in [0.02, 50]$ (identifying the thick black edge separating white "no patterns" left from shaded "pattern-forming region" right).
  - *Calibration (1200 dpi):* Log10 x-axis borders + major tick ($10^0, 10^1, 10^2$, residual `0.33 px`); log10 y-axis borders + 3 major ticks ($10^{-2}..10^2$, residual `0.20 px`).
  - *Reproduction verdicts:* Swept directly on existing model `single_species_stability_boundary_hill`: `cmp_fig3b_stability_boundary` **reproduced** (11/11 points reproduced; nose at $c_0 = 1$ reads `4.001` vs. exact `4.000`).

- **Fig. 4 (p. 3) — Left out:**
  - *Reason:* 3D oblique perspective waterfall surface plot of transient PDE dynamics $c(x/\ell, t/\tau_D)$ for $L/\ell = 3\pi$, with overlapping occluded slices and perspective projection lacking orthogonal calibrated 2D axes. Its stated quantitative properties (fastest-growing mode $n = 2$ and final single-peak steady state) are already checked in `cmp_fastest_mode_fig4` and `cmp_final_peak_count_fig4`.

- **Fig. 5(a) (p. 4) — Read (`fig5a_boundary_I_II`):**
  - *Curves:* ASDM linear stability diagram at $\mathrm{Pe} = 0$, $\rho_a = 1$: Region I/II Turing boundary $D_a/D_s$ read horizontally at 8 values of $\rho_s \in [1.05, 1.95]$ (avoiding $\rho_s = 1.5$ where the `+` marker touches the boundary).
  - *Calibration (1200 dpi):* Borders + central x-tick ($D_a/D_s = 0, 1, 2$, residual `0.17 px`); borders + $\rho_s = 1$ horizontal line ($\rho_s = 0, 1, 2$, residual `0.00 px`).
  - *Additional calibrated measurements:* The unstated `+` marker centroid measures $(D_a/D_s, \rho_s) = (0.0752, 1.5000)$ (settling `q_fig5_marker_coords` via `q_fig5_markers_calibrated`), and the horizontal Region III/V boundary measures $\rho_s = 0.1717$ (matching $3 - 2\sqrt{2} = 0.1716$). Lower-$\rho_s$ boundaries (Regions III/IV/V) are left out of `figures.yaml` because existing `asdm_stability_boundary` covers only $\rho_s > \rho_a$.
  - *Reproduction verdicts:* Swept directly on existing model `asdm_stability_boundary`: `cmp_fig5a_boundary_I_II` **reproduced** (8/8 points reproduced).

- **Fig. 5(b) (p. 4) — Read (`fig5b_boundary_I_II`):**
  - *Curves:* ASDM linear stability diagram at $\mathrm{Pe} = 7$, $\rho_a = 1$: Region I/II active-fluid stability boundary $D_a/D_s$ read horizontally at 10 values of $\rho_s \in [1.05, 1.95]$.
  - *Calibration (1200 dpi):* Borders + central x-tick ($D_a/D_s = 0, 1, 2$, residual `0.17 px`); borders + left y-tick ($\rho_s = 0, 1, 2$, residual `0.00 px`).
  - *Additional calibrated measurements & finding (`q_fig5b_finite_k_grid`, `q_fig5_markers_calibrated`):* The `×` marker centroid measures $(D_a/D_s, \rho_s) = (0.5002, 1.5002)$, confirming $(0.5, 1.5)$. However, the vertical Region I/III boundary lies at $D_a/D_s = 0.7082$ (instead of the continuum small-$k$ limit $\mathrm{Pe}/4 - 1 = 0.75$) and the bottom boundary of Region III slopes from $\rho_s = 0.2220$ ($D_a/D_s = 0.8$) down to $0.2097$ ($D_a/D_s = 1.9$) instead of lying flat at $3 - 2\sqrt{2} = 0.1716$. Both match Eq. (15) on a discrete $k$-grid with minimum nonzero wavenumber $k_1\ell = \pi/20 \approx 0.1571$ (predicting $D_a/D_s = 0.7079$ and $\rho_s = 0.2226 \to 0.2101$).
  - *Reproduction verdicts:* Swept directly on existing model `asdm_stability_boundary`: `cmp_fig5b_boundary_I_II` **reproduced** (10/10 points reproduced).

- **Fig. 5(c) (p. 4) — Read (`fig5c_Pe50`, `fig5c_Pe25`, `fig5c_Pe0_cross`, `fig5c_Pe0_plus`):**
  - *Curves:* ASDM dispersion relations $\tau_D \lambda(k)$ vs. $k\ell$ ($\tau_D = \ell^2/D_s$): solid curves at `×` parameters $(0.5, 1.5)$ for $\mathrm{Pe} = 50$ (11 points), $\mathrm{Pe} = 25$ (9 points), and $\mathrm{Pe} = 0$ (8 points); dashed reference curve at `+` parameters $(0.075, 1.5)$ for $\mathrm{Pe} = 0$ (9 points). Thin zero-line and serif curve labels ("Pe=50", "Pe=25", "Pe=0") are removed via morphological opening and band filtering.
  - *Calibration (1200 dpi):* Borders + 2 interior x-ticks ($k\ell = 0, 2, 4, 6$, residual `0.00 px`); borders + 2 interior y-ticks ($\tau_D\lambda = -5, 0, 5, 10$, residual `0.20 px`).
  - *Reproduction verdicts:* Swept directly on existing models `asdm_dispersion_eq15` and `asdm_dispersion_pe0`: `cmp_fig5c_Pe50` **reproduced**, `cmp_fig5c_Pe25` **reproduced**, `cmp_fig5c_Pe0_cross` **reproduced**, `cmp_fig5c_Pe0_plus` **reproduced** (4/4 reproduced; e.g. $\mathrm{Pe}=50$ peak at $k\ell=2$ reads `9.464` vs. `9.468`, confirming $\rho_s = 1.5$ and $D_a/D_s = 0.5$ at `×` and $0.075$ at `+`).

- **Fig. 5(d) (p. 4) — Read (`fig5d_a`, `fig5d_s`, `fig5d_v`):**
  - *Curves:* ASDM nonhomogeneous 3-peak steady state on $x/\ell \in [0, 20]$ with periodic BCs: activator $a(x/\ell)$ (solid black, left axis, 14 points), substrate $s(x/\ell)$ (solid gray, left axis, 14 points), and flow velocity $v/U(x/\ell)$ (dashed black, right axis, 15 points). The gray substrate curve (`gray ~ 149`) is isolated by intensity thresholding plus `5x5` disk opening to eliminate black-stroke antialiasing fringes.
  - *Calibration (1200 dpi):* Borders + 3 interior x-ticks ($x/\ell = 0, 5, 10, 15, 20$, residual `0.20 px`); 5 left-border y-ticks ($a, s = 1..5$, residual `0.20 px`); 3 right-border y-ticks ($v/U = -0.2, 0, 0.2$, residual `0.00 px`).
  - *Reproduction verdicts & major finding (`q_fig5d_substrate_rho_s_mismatch`):* Evaluated via curve-port wrapper `asdm_pde_profile_eq14` (over existing `asdm_pde_eq14`, using `n_grid = 128`, `dt = 0.004`, centered at $x/\ell = 10$) at the stated `×` parameters ($D_a/D_s = 0.5$, $\rho_a = 1$, $\rho_s = 1.5$, $\mathrm{Pe} = 7$):
    - `cmp_fig5d_v`: **reproduced** (15/15 points reproduced).
    - `cmp_fig5d_a`: **reproduced** (peak prediction `5.166` lies slightly below raw reading `5.368 [5.200, 5.536]`, falling inside the printed-input rounding band `[4.635, 5.698]` for integer `Pe = 7`, `rho_a = 1`).
    - `cmp_fig5d_s`: **not reproduced** (10/14 points not reproduced: Eq. 14 at $\rho_s = 1.5$ predicts substrate plateau $s_{\max} = 2.038$ [rounding band `2.007–2.069`], whereas Fig. 5(d) plots $s_{\max} = 1.400 \pm 0.052$ and $s_{\min} = 0.126$ vs. $0.100$ — a 46% mismatch). Diagnostic tests show that changing only $\rho_s$ from `1.5` to `1.0` matches all three plotted curves simultaneously to `< 0.5%` ($s_{\max} = 1.402$, $s_{\min} = 0.125$, $a_{\max} = 5.424$, peak half-width at $a = 5$ is `0.197` vs. measured `0.201`, $(v/U)_{\max} = 0.253$), strongly suggesting Fig. 5(d) was generated with $\rho_s = 1$ rather than the `(0.5, 1.5)` `×` point used in Fig. 5(b)–(c).
