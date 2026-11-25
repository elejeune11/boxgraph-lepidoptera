# Figure Digitization Log: `tjhung2018`

Calibrated figure digitization of Tjhung, Nardini, and Cates, *Phys. Rev. X* **8**, 031080 (2018) (PDF pp. 1–19, `031080-1` to `031080-19`, including Appendices A–G).

- **Total calibrated figure readings**: **32** (all `mode: code` in `figures.yaml`, self-contained fast readers in `figures/*.py`, saved in `readings.json`).
- **Total digitized data points**: **345** (`0` missing/errored points).
- **Visual verification**: Every single one of the 32 `boxgraph figure check` overlays was rendered under `figwork/<id>.png`, compiled into 6 high-visibility contact sheets (`figwork/sheet_vis_1.png`–`sheet_vis_6.png`), and visually inspected to confirm 100% alignment of every magenta mark with its intended curve, marker, contour, or arrow.
- **New models & sweep comparisons**: Added 3 Fig. 2 parameterization models (`binodals_flat_combo_eq9`, `binodals_flat_combo_eqB2`, `binodals_flat_combo_eqB2_corrected`), 1 reconstructed 2D log-screened two-droplet ODE (`two_droplet_ode_2d_log_screened`), 27 sweep runs, and 32 curve-level comparisons testing printed, `variant_of`, and reconstructed models directly against the digitized curves.

---

## Panel-by-Panel Inventory & Verification

### Fig. 1 (PDF p. 2, `031080-2`)
- **Fig. 1(a) & Fig. 1(b)** — **Left out**: Qualitative 2D scalar-field color snapshots $\phi(\mathbf{r}, t)$ at four time points each (`t = 38500, 39700, 40500, 40950` for bubbly phase separation at $\zeta=4, \lambda=1, \phi_0=-0.4$; `t = 100, 400, 40000, 60000` for arrested cluster phase at $\zeta=-4, \lambda=-1, \phi_0=-0.6$). No coordinate axes, curves, points, or quantitative data plots to digitize.

### Fig. 2 (PDF p. 6, `031080-6`)
- **Fig. 2 (flat-interface binodals $\phi_1, \phi_2$ vs $\zeta - 2\lambda$)** — **Read**:
  - **Readings**: `fig2_phi1` (14 pts, lower magenta curve $\phi_1$, `fig2_phi1_curve`) and `fig2_phi2` (14 pts, upper blue curve $\phi_2$, `fig2_phi2_curve`), sampled at $\zeta - 2\lambda \in \{-5.8, -5, -4, -3, -2, -1, -0.2, 0.2, 1, 2, 3, 4, 5, 5.8\}$ (`1200 dpi`, crop `[1040, 740, 4800, 3580]`).
  - **Calibration**: Dynamic inward-tick detection along top edge (7 x-ticks `-6..6`) and right edge (7 y-ticks `-1.5..1.5`); max calibration residual **0.43 px** on both x and y (`< 0.001` in $\phi$).
  - **Verdicts & settled questions**:
    - `binodals_flat_combo_eq9` and `binodals_flat_combo_eqB2_corrected` (at $\lambda = 1$) both achieve **14/14 `reproduced`** on both $\phi_1$ and $\phi_2$ with raw discrepancies `< 0.0008` (e.g., at $\zeta - 2\lambda = 2$: digitized `(-1.095, 0.8634)` vs exact `(-1.0956, 0.8637)`).
    - Printed `binodals_flat_combo_eqB2` (`-6*zeta*K^2*phi`, evaluated at $\lambda = 1$) fails catastrophically (`not reproduced` on both curves even after automatic `input_rounding`; 14/14 raw `not reproduced` on $\phi_1$ and 12/14 raw `not reproduced` + 2 `close` on $\phi_2$), quantitatively settling `q_eqB2_misprint_6zeta` across the entire plotted range.
  - **Note on dashed spinodals**: Horizontal reference lines $\phi_s = \pm 1/\sqrt{3} \approx \pm 0.577$ are analytical constants already represented in `model_b_bulk_eq3` (`phi_s_pos_eq3`) and left as reference geometry.

### Fig. 3 (PDF p. 6, `031080-6`)
- **Fig. 3(a) (outer curved-droplet density $\phi_-(R)$ vs $R$)** — **Read**:
  - **Readings**:
    - `fig3a_theory_A` (11 pts, red solid analytical curve, $\zeta=-1, \lambda=0.5$, $R \in [3.5, 38]$, `fig3a_phi_minus_theory_A`).
    - `fig3a_theory_B` (9 pts, blue solid analytical curve, $\zeta=-4, \lambda=-1$, $R \in [5, 38]$, `fig3a_phi_minus_theory_B`).
    - `fig3a_sim_A` (5 pts, red `+` simulation markers at $R \in [8.90, 34.56]$, `fig3a_phi_minus_sim_A`, `origin: run`), detected via 2× downsampled + full-resolution `+` matched-filter correlation.
    - `fig3a_sim_B` (10 pts, blue `×` simulation markers at $R \in [6.79, 34.37]$, `fig3a_phi_minus_sim_B`, `origin: run`), cleanly separating both closely spaced touching pairs `(8.25, 8.74)` and `(10.42, 11.06)` via `×` matched-filter correlation.
  - **Calibration**: `1200 dpi`, crop `[5380, 740, 9280, 3460]`; 9 x-ticks (`0..40` step 5) and 4 y-ticks (`-0.7..-1.0` step `0.1`), residual **0.35 px** (x) and **0.25 px** (y). Inset and legend regions masked out explicitly.
  - **Inset**: Qualitative 2D droplet density heatmap showing $\phi_+$, $\phi_-$, and $2R$ — left out as schematic illustration.
- **Fig. 3(b) (inner curved-droplet density $\phi_+(R)$ vs $R$)** — **Read**:
  - **Readings**:
    - `fig3b_theory_A` (10 pts, red solid analytical curve, $R \in [4.5, 38]$, `fig3b_phi_plus_theory_A`).
    - `fig3b_theory_B` (7 pts, blue solid analytical curve, $R \in [13.5, 38]$, `fig3b_phi_plus_theory_B`).
    - `fig3b_sim_A` (5 pts, red `+` simulation markers at $R \in [8.90, 34.57]$, matching Fig. 3a radii to `< 0.01`, `fig3b_phi_plus_sim_A`, `origin: run`).
    - `fig3b_sim_B` (3 pts, blue `×` simulation markers at $R = 14.10, 26.04, 34.38$, matching Fig. 3a radii to `< 0.005`, `fig3b_phi_plus_sim_B`, `origin: run`).
  - **Calibration**: `1200 dpi`, crop `[5380, 3480, 9280, 6200]`; 9 x-ticks (`0..40`) and 4 y-ticks (`1.3..1.0`), residual **0.36 px** (x) and **0.58 px** (y).
- **Verdicts & settled questions (Fig. 3a–b)**:
  - Directly settles `q_appC_eqC2_signs_and_limits` and `q_audit_appC_scope_and_w2_clipping`:
    - Printed Eq. (C2) (`curved_binodals_appC`) fails catastrophically on all 4 analytical curves in raw comparison (errors up to `-0.074` on Fig. 3a A at $R=5$, non-monotonic `-0.15` on Fig. 3b A at $R=5$, and `+0.091` on Fig. 3b B at $R=16$); `cmp_fig3b_appC_printed_A` remains **`not reproduced`** even after BoxGraph's wide ±50% automatic `input_rounding` bands on integer parameters (`K=1, d=2, zeta=-1`).
    - Corrected Eq. (C2) (`curved_binodals_appC_corrected`) matches Fig. 3a A & B and Fig. 3b A to `< 0.001` for $R \ge 12$, and improves Fig. 3b B by `> 10×` (residual `0.001–0.009`), confirming that the plotted Fig. 3b B theory curve sits tightly between 1-step (`1.2484` at $R=16$) and 2-step/fixed-point (`1.2617`/`1.2664` vs digitized `1.2571`) Appendix C iteration.
    - Leading $\mathcal{O}(1/R)$ (`curved_binodals_O1R_eq16_19`) underpredicts small-$R$ curvature shifts in Region B by `0.033–0.035`.
    - Comparing `curved_binodals_appC_corrected` against the 23 digitized simulation markers confirms 3/5 exact matches on Fig. 3b A ($R \ge 15.6$) and systematic `0.004–0.015` finite-width/finite-box offsets at smaller $R$ and in Fig. 3a.

### Fig. 4 (PDF p. 8, `031080-8`)
- **Fig. 4(a) (mean-field phase diagram in $(\zeta, \lambda)$)** — **Read**:
  - **Readings**:
    - `fig4a_boundary_AB` (11 pts, black $\sigma = 0$ curve separating Regions A & B, $\zeta \in [-3.9, -0.4]$, `fig4a_boundary_AB_curve`).
    - `fig4a_boundary_AC` (11 pts, magenta $\sigma_{\text{bubble}} = 0$ curve separating Regions A & C, $\zeta \in [0.4, 3.9]$, `fig4a_boundary_AC_curve`).
    - `fig4a_sim_forward` (6 pts, blue circle forward-Ostwald markers at `(0.20, 1.10), (0.80, 0.40), (1.00, -0.50), (1.40, -0.30), (1.80, -0.10), (2.00, 0.00)`, `fig4a_sim_forward_pts`, `origin: run`).
    - `fig4a_sim_reverse` (5 pts, red square reverse-Ostwald markers at `(0.80, 1.40), (1.30, 0.65), (1.60, 0.80), (2.40, 0.20), (4.00, 1.00)`, handling the half-clipped border square at $\zeta=4.0$, `fig4a_sim_reverse_pts`, `origin: run`).
  - **Calibration**: `1200 dpi`, crop `[880, 740, 4980, 4560]`; 9 x-ticks (`-4..4`) and 9 y-ticks (`-2..2`), residual **0.39 px** (x) and **0.25 px** (y). Three heatmap inset boxes excluded from masks.
  - **Verdicts**: `phase_boundary_eq23` achieves **11/11 `reproduced`** on `cmp_fig4a_boundary_AB` and **11/11 `reproduced`** on `cmp_fig4a_boundary_AC` without input-rounding help (max raw discrepancy `~0.002–0.006` on the steepest segments). Independent check of all 11 simulation markers against `pseudotension_flat_eq20` confirms 6/6 blue circles have $\sigma_{\text{bubble}} > 0$ (`+0.047` to `+0.260`) and 5/5 red squares have $\sigma_{\text{bubble}} < 0$ (`-0.034` to `-0.357`).
  - **Insets**: Three schematic 2D two-droplet/two-bubble snapshots — left out (no data axes).
- **Fig. 4(b) (two-droplet radii $R_1(t), R_2(t)$ vs $\log_{10} t$)** — **Read**:
  - **Readings (all 8 curves, 90 pts total)**:
    - Region A ($\zeta=-1, \lambda=0.5$, red): `fig4b_theory_A_R1` (11 pts, dashed, $t \in [200, 16500]$), `fig4b_theory_A_R2` (10 pts, dashed, $t \in [200, 16500]$), `fig4b_sim_A_R1` (11 pts, solid, $t \in [200, 50000]$), `fig4b_sim_A_R2` (11 pts, solid, $t \in [200, 23500]$).
    - Region B ($\zeta=-4, \lambda=-1$, blue): `fig4b_theory_B_R1` (12 pts, dashed, $t \in [200, 150000]$), `fig4b_theory_B_R2` (12 pts, dashed, $t \in [200, 150000]$), `fig4b_sim_B_R1` (12 pts, solid, $t \in [200, 100000]$), `fig4b_sim_B_R2` (12 pts, solid, $t \in [200, 100000]$).
  - **Calibration**: `1200 dpi`, crop `[960, 4520, 4980, 7340]`; dynamic detection of the 4 major decade ticks `100, 1000, 10000, 100000` (residual **0.35 px**) revealed that the right plot border is **$t = 200{,}000$** (not `100000`); linear y-axis calibrated from 6 ticks `0..25` (residual **0.34 px**). Dashed curves interpolated across dash gaps via local linear fits after separating the continuous solid band.
  - **Verdicts & settled questions**: Directly settles `q_appE_mass_conservation_d2_vs_d3` and `q_audit_fig4b_2d_log_factor_and_fig6_scaling`:
    - Verbatim printed Eqs. (E1)–(E2) (`two_droplet_ode_eqE1_E2`) fails on 10/10 raw points of Region A $R_2$ (evaporating $R_2$ by $t \approx 10{,}000$ vs `~17500` in Fig. 4b red dashed, ~1.7× too fast), caps Region A $R_1 \le 19.07$ by area conservation (vs `> 20.0` in Fig. 4b red dashed), and converges Region B to `13.46` (vs `13.26` in Fig. 4b blue dashed).
    - Reconstructed 2D log-screened ODE (`two_droplet_ode_2d_log_screened`, dividing $\dot{R}_i$ by $\ln(\ell/R_i)$ with $\ell = 50$) simultaneously resolves all three raw discrepancies — matching the ~1.7× slowdown, the Region A $R_1 > 19.1$ overshoot (`20.1` at $t = 16500$), and the Region B $R_{\text{eq}}$ undershoot (`13.1–13.2` for $\ell = 50$–$64$).
  - **Inset**: Schematic 2D two-droplet snapshot indicating $2R_1(t), 2R_2(t)$ — left out.

### Fig. 5 (PDF p. 10, `031080-10`)
- **Fig. 5(a) (steady-state phase diagram vs global density $\phi_0$)** — **Left out as image reading (printed labels already recorded)**: 2×4 matrix of qualitative 2D stochastic simulation snapshots at $\phi_0 = -1.1, -0.4, 0.6, 1.0$ with a schematic bottom $\phi_0$ arrow. Its exact printed binodal labels $\phi_1 = -1.096$ and $\phi_2 = 0.864$ are already recorded in `data.yaml` (`phi1_z2l2_fig5a`, `phi2_z2l2_fig5a`) and compared against `binodals_flat_eq9`, `binodals_flat_eqB2`, `binodals_flat_eqB2_corrected`, and `amb_plus_1d_flat_pde`.
- **Fig. 5(b) (bubble radius probability distribution $P(R)$ vs $R$)** — **Read**:
  - **Readings (4 curves, 57 pts total, `origin: run`)**:
    - `fig5b_PR_t512` (13 pts, red, $R \in [2.75, 5.5]$, peak `0.224` at $R = 3.85$).
    - `fig5b_PR_t4096` (14 pts, green, $R \in [3.5, 7.0]$, peak `0.204` at $R = 4.35$).
    - `fig5b_PR_t32768` (15 pts, blue, $R \in [4.0, 8.0]$, peak `0.150` at $R = 5.25$).
    - `fig5b_PR_t131072` (15 pts, magenta, $R \in [4.0, 8.0]$, peak `0.148` at $R = 5.5$, confirming steady state reached by $t \simeq 32000$).
  - **Calibration**: `1800 dpi`, crop `[9870, 1110, 13350, 4050]`; 10 x-ticks (`1..10`) and 6 y-ticks (`0..0.25`), residual **0.44 px** (x) and **0.40 px** (y). Legend excluded.
  - **Finding**: All four curves integrate to $\int P(R)\,dR = 0.250$, proving discrete histogram binning with $\Delta R = 0.25$ and $\sum P(R_i) = 1$; steady-state first moment $\bar{R}(t=32768) = 5.56$ matches Fig. 5(c) at $D = 0.200$ ($\bar{R} = 5.663$).
- **Fig. 5(c) (mean bubble radius $\bar{R}$ and bubble number density $\bar{n}$ vs noise $D$)** — **Read**:
  - **Readings (2 dual-axis marker sets, 22 pts total, `origin: run`)**:
    - `fig5c_Rbar` (11 red `+` markers, left y-axis $\bar{R} \in [2.34, 12.48]$ across $D \in \{0.12, 0.13, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.80\}$).
    - `fig5c_nbar` (11 blue `+` markers, right y-axis $\bar{n} \in [0.000239, 0.00505]$ across the exact same 11 noise strengths $D$, cleanly separating the overlapping blue pair at $D = 0.12, 0.13$ via matched-filter correlation).
  - **Calibration**: `1800 dpi`, crop `[9930, 3990, 13920, 6990]`; 6 x-ticks (`0..1.0`, residual **0.40 px**), 8 left y-ticks (`0..14` — noting top frame border sits at $\bar{R} = 15$, residual **0.36 px**), and 9 right y-ticks (`0..0.008`, residual **0.44 px**).

### Fig. 6 (PDF p. 16, `031080-16`)
- **Fig. 6 (two-droplet interface contours and vapor current field $\mathbf{J}(x, y)$ at $t = 1000$)** — **Read**:
  - **Readings (4 calibrated datasets, 54 pts total, `origin: run`)**:
    - `fig6_droplet_contours_normal` (16 pts, top half $y > 30$, Normal Ostwald $\zeta=-0.8, \lambda=-0.4$, blue $\phi = 0$ contour $y(x)$ across both droplets: circle fits give $R_1 = 23.25$ at $x_1 = 32.10$ and $R_2 = 5.85$ at $x_2 = 95.03$).
    - `fig6_droplet_contours_reverse` (16 pts, bottom half $y < 30$, Reverse Ostwald $\zeta=-1.6, \lambda=-0.8$, blue $\phi = 0$ contour $y(x)$: circle fits give $R_1 = 23.25$ at $x_1 = 32.10$ and $R_2 = 7.53$ at $x_2 = 95.03$).
    - `fig6_Jx_profile_normal` (11 pts, top-half midline vapor arrows along $y \approx 32$, normalized signed $J_x(x)/|J_x(68)|$ showing leftward flow `-1.28..-0.97` in the inter-droplet gap $x \in [59.6, 83.5]$ and rightward periodic flow `+0.99..+1.04` in outer vapor $x \in [112.4, 124.4]$).
    - `fig6_Jx_profile_reverse` (11 pts, bottom-half midline vapor arrows along $y \approx 28$, normalized signed $J_x(x)/|J_x(68)|$ showing rightward flow `+0.97..+1.52` in the inter-droplet gap $x \in [57.2, 80.4]$ and leftward periodic flow `-1.10..-0.91` in outer vapor $x \in [111.7, 123.8]$).
  - **Calibration**: `1500 dpi`, crop `[6950, 925, 11325, 3275]`; 7 x-ticks (`0..120`, residual **0.37 px**) and 7 y-ticks (`0..60`, residual **0.20 px**), revealing the actual plotted lattice extent is **$x \in [0, 128]$, $y \in [0, 64]$** ($128 \times 128$ periodic simulation domain).
  - **Verdicts & settled questions**: Directly settles `q_fig6_vs_fig4b_parameter_mismatch` and `q_unstated_box_sizes_and_radii`, and an independent $128 \times 128$ 2D AMB+ PDE simulation from $R_1(0) = 23.25, R_2(0) = 6.7$ reproduces both $R_2(1000)$ (`5.56` normal vs `7.51` reverse), $|J_x(68)| \approx 3.1$–$3.7 \times 10^{-4}$ (matching `"around 10^-4"`), and the exact sign reversal across both vapor regions.

### Appendices A–G (PDF pp. 14–17)
- Only **Fig. 6** (p. 16, Appendix F) appears across Appendices A–G; no other plots exist in the paper.
