# Figure Readings — `haviv2008`

Calibrated `mode: code` digitization of every quantitative plot across all 10 PDF pages of Haviv, Gov, Ideses & Bernheim-Groswasser, *Eur. Biophys. J.* **37**, 447–454 (2008; DOI `10.1007/s00249-007-0236-1`), `corpus/actin_bundling/2008_Haviv.pdf` (pp. 1–8 main article, pp. 9–10 Supporting Information converted from `249_2007_236_MOESM1_ESM.doc`). All four quantitative panels (**Fig. 2a–d**, p. 4 / 450) were rendered at **1200 dpi** (plot boxes 2408–2568 × 2488–2536 px; native embedded bitmap resolution ~150 dpi, i.e., 8 px per native pixel, so ±4.0 px native quantization is included in every reading interval). Axes were calibrated purely from detected neutral-black frame spines and inward tick marks via `figures/_calib.py` (residuals `0.13–3.54 px` across `~2500 px`, `< 0.14%`), and all **14 overlay PNGs** in `figwork/` were generated with `boxgraph figure check` and visually inspected via compact contact sheets to verify every magenta cross and vertical error-bar span.

- **Pages 1–2 (447–448) — No figures:** Title, abstract, introduction, materials & methods, and onset of results (text only).

- **Fig. 1a–d (p. 3 / 449) — Left out:**
  - *Reason:* High-resolution cryo-TEM micrographs of actin/fascin systems at $[A] = 5\,\mu\text{M}$ and $[F]/[A] = 1/200$ (a, individual filaments), $1/2$ (b), $1/1$ (c), and $6/1$ (d) (scale bars $0.2\,\mu\text{m}$); micrographs with no plotted numeric axes or data curves.

- **Fig. 1e (p. 3 / 449) — Left out:**
  - *Reason:* Schematic illustration of cylindrical bundle cross-section showing bundle diameter $D_b$ and effective filament center-to-center distance $D_f$, plus side cartoon of yellow actin filaments cross-linked by blue fascin spheres; drawing with no numeric axes.

- **Fig. 2a (p. 4 / 450) — Read (`fig2a_data`, `fig2a_fit`):**
  - *Series:* Linear plot of normalized bundle-size distribution $P(N_b)$ vs. filament number $N_b$ at $[F]/[A] = 3/1$, $[A] = 5\,\mu\text{M}$:
    1. `fig2a_data` (`P_fig2a_measured`, `origin: measured`): all **29** evenly spaced histogram bins (blue filled circles, spacing $\Delta N_b \approx 15.6$ from $N_b = 11.0$ to $448.2$, splitting the adjacent zero-count double marker at $N_b = 385.9, 400.5$): **19 nonzero bins** (`0.781, 0.999, 0.812, 0.906, 0.655, 0.532, 0.469, 0.220, 0.218, 0.220, 0.094, 0.063, 0.066, 0.030, 0.033, 0.063, 0.033, 0.032, 0.032`) and **10 zero-count bins** sitting on the bottom $x$-axis ($N_b = 245.3, 276.2, 307.6, 323.1, 338.8, 369.6, 385.9, 400.5, 416.9, 433.2$).
    2. `fig2a_fit` (`P_fig2a_eq1_fit`, `origin: predicted`): peak-normalized Eq. (1) solid black fit curve sampled at **13** unobstructed $N_b$ positions (`15, 35, 50, 65, 80, 95, 110, 130, 160, 175, 190, 237, 250`), spanning $P = 0.995$ near the peak down to $0.019$ at $N_b = 250$.
  - *Calibration (1200 dpi):* Linear $x$-axis from left spine + 4 clean inward ticks (`0, 50, 100, 150, 200`, residual `0.78 px`, right spine at `N_b = 448.2`); linear $y$-axis from bottom spine + 4 clean inward ticks (`0.0, 0.2, 0.4, 0.8, 1.0`, residual `0.83 px`; tick `0.6` occluded by rising curve).
  - *Verdicts:*
    - `cmp_eq1_vs_fig2a_fit` (`reproduction`, `size_distribution_eq1` with Table 1 $\phi = 3/1$ row $2k_{\text{on}}/D = 0.4396$, $k_{\text{off}}/D = 0.04$, $f = 1.0$): **reproduced (13/13 points)**, confirming peak-height normalization $R_n = f \exp(-(k_{\text{on}}/D)^2/(k_{\text{off}}/D))$ (`q_unstated_normalization_Rn_and_nonzero_P0`).
    - `cmp_eq1_vs_fig2a_data` (`fit_quality`): **not reproduced (25/29 reproduced, 4/29 not reproduced)**; the 4 departures are the noisy heavy-tail nonzero bins at $N_b = 260.3$ (`0.063` vs. `0.011`), $292.2$ (`0.033` vs. `0.0046`), $354.3$ (`0.032` vs. `0.0008`), and $448.2$ (`0.032` vs. `0.0001`).

- **Fig. 2b (p. 4 / 450) — Read (`fig2b_data`, `fig2b_fit`):**
  - *Series:* Semi-log plot of `Log(P(N_b))` vs. $N_b$ at $[F]/[A] = 3/1$. The ordinate is **natural logarithm** $\ln P(N_b)$ (not $\log_{10}$); the reader scripts return linear $P = \exp(\ln P)$ paired with `scale: ln` data in `data.yaml` so BoxGraph compares directly in natural-log space:
    1. `fig2b_data` (`lnP_fig2b_measured`, `origin: measured`, `scale: ln`): all **19** isolated nonzero blue circles ($N_b = 11.5\text{–}448.1$, $\ln P = 0.002\text{ to }-3.468$). Every single bin matches $\ln(k/32)$ to `0.002–0.010` for exact integer counts $k = [25, 32, 26, 29, 21, 17, 15, 7, 7, 7, 3, 2, 2, 1, 1, 2, 1, 1, 1]$, proving **exact sample size $n = \sum k = 200$ bundles** at $\phi = 3/1$! Independent cross-check: the binned sample moments of these 200 bundles give $N_{b,\text{mean}} = 73.41 \pm 4.38$ and $D_{b,\text{mean}} = 71.49 \pm 2.05\,\text{nm}$, agreeing with Supporting Table S1 ($74.0 \pm 4.4$, $71.8 \pm 2.1\,\text{nm}$) well inside 1 SE.
    2. `fig2b_fit` (`lnP_fig2b_eq1_fit`, `origin: predicted`, `scale: ln`): Eq. (1) solid black curve sampled at **15** unobstructed $N_b$ positions (`15, 35, 50, 65, 80, 95, 112, 128, 145, 175, 205, 225, 245, 270, 280`), spanning $\ln P = -0.010$ at $N_b = 35$ down to $-5.071$ at $N_b = 280$.
  - *Calibration (1200 dpi):* Linear $x$-axis from left spine + 8 inward ticks (`0, 50, ..., 400`, residual `1.11 px`, right spine at `N_b = 448.2`); linear $y$-axis in $\ln P$ from 10 inward ticks (`0, -1, -1.5, ..., -5`, residual `1.53 px`; tick `-0.5` occluded by rising curve).
  - *Verdicts:*
    - `cmp_eq1_vs_fig2b_fit` (`reproduction`): **reproduced (15/15 points)** in natural-log space, decisively confirming $\ln P(N_b)$ (`q_unstated_normalization_Rn_and_nonzero_P0`).
    - `cmp_eq1_vs_fig2b_data` (`fit_quality`): **not reproduced (9/19 reproduced, 10/19 not reproduced)**; matches closely across the peak and upper flank ($N_b \le 136$), then diverges systematically in the finite-sample heavy tail ($N_b = 261, 292, 354, 448$: $\ln P = -2.77\text{ to }-3.47$ vs. Eq. (1) $-4.55, -5.38, -7.11, -9.83$), quantifying the caption's remark that large-width statistics are noisy.

- **Fig. 2c (p. 4 / 450) — Read (`fig2c_koff_data`, `fig2c_kon_data`, `fig2c_koff_arrhenius_curve`, `fig2c_kon_constant_curve`, `fig2c_kon_adsorption_curve`, `fig2c_kon_coat_hindered_curve`):**
  - *Double mislabeling confirmed (`q_fig2c_axis_and_caption_mislabeling`, `q_fig2_digitization_settles_prefactors_ln_mislabeling_and_bars`):* The caption swaps both color and symbol ("black squares and red circles") and the $y$-axis prints `k_on/D, 10k_off/D` omitting the factor of 2. Calibrated digitization proves:
    1. `fig2c_koff_data` (`koff_over_D_fig2c_data`, `origin: fitted`, `error_bars: range`): **5** black hollow circles with black vertical error bars (axis plots $10 k_{\text{off}}/D$; reader divides by 10 to report $k_{\text{off}}/D$):
       - $\phi = 0.223$: $0.09014$ (Table 1: `0.09`, `+0.16%`), bar `[0.07214, 0.1082]` (`-20.0% / +20.1%`)
       - $\phi = 0.507$: $0.08001$ (Table 1: `0.08`, `+0.02%`), bar `[0.06403, 0.09596]` (`-20.0% / +19.9%`)
       - $\phi = 0.996$: $0.06207$ (Table 1: `0.062`, `+0.11%`), bar `[0.04964, 0.07438]` (`-20.0% / +19.8%`)
       - $\phi = 2.996$: $0.04011$ (Table 1: `0.04`, `+0.29%`), bar `[0.03210, 0.04807]` (`-20.0% / +19.8%`)
       - $\phi = 6.009$: $0.01909$ (Table 1: `0.019`, `+0.48%`), bar `[0.01473, 0.02280]` (`-22.8% / +19.5%`)
    2. `fig2c_kon_data` (`two_kon_over_D_fig2c_data`, `origin: fitted`, `error_bars: range`): **5** orange-filled red squares with red vertical error bars plotting **$2 k_{\text{on}}/D$** (NOT $k_{\text{on}}/D$):
       - $\phi = 0.203$: $0.4917$ (Table 1: `0.4962`, `-0.90%`), bar `[0.3947, 0.5964]` (`-19.7% / +21.3%`)
       - $\phi = 0.500$: $0.5723$ (Table 1: `0.5724`, `-0.01%`), bar `[0.4561, 0.6894]` (`-20.3% / +20.5%`)
       - $\phi = 0.988$: $0.5724$ (Table 1: `0.5695`, `+0.50%`), bar `[0.4561, 0.6841]` (`-20.3% / +19.5%`)
       - $\phi = 3.004$: $0.4434$ (Table 1: `0.4396`, `+0.86%`), bar `[0.3543, 0.5263]` (`-20.1% / +18.7%`)
       - $\phi = 5.997$: $0.2159$ (Table 1: `0.2159`, `-0.01%`), bar `[0.1719, 0.2596]` (`-20.4% / +20.3%`)
    3. `fig2c_koff_arrhenius_curve` (`koff_over_D_fig2c_arrhenius`, `origin: predicted`): solid black Arrhenius curve sampled at **13** $\phi$ values (`0.05..9.80`, divided by 10; red dotted curve coalesces along its upper edge for $\phi \ge 3$).
    4. `fig2c_kon_constant_curve` (`two_kon_over_D_fig2c_constant`, `origin: predicted`): solid red horizontal constant on-rate line sampled at **8** $\phi$ values (`0.5..9.5`), sitting flat at $2k_{\text{on}}/D = 0.4561$ (essentially the arithmetic mean `0.4587` of Table 1's 5 rows).
    5. `fig2c_kon_adsorption_curve` (`two_kon_over_D_fig2c_adsorption`, `origin: predicted`): dashed red Langmuir adsorption curve sampled on **15** isolated dashes (`0.36..9.60`).
    6. `fig2c_kon_coat_hindered_curve` (`two_kon_over_D_fig2c_coat_hindered`, `origin: predicted`): dotted red coat-hindered curve sampled at **12** $\phi$ values (`0.35..9.50`; coalesced with the solid black curve for $\phi \ge 3$ because $1.054\,p_{fa}(\phi)\approx 1.00$ there).
  - *Calibration (1200 dpi):* Linear $x$-axis from left/right spines + 9 inward ticks (`0..10`, 11 marks, residual `0.91 px`); linear $y$-axis from bottom spine + 4 clean neutral inward ticks (`0.0, 0.4, 0.6, 0.8, 1.0`, residual `0.13 px`; tick `0.2` occluded by rising red curves).
  - *Verdicts & unstated prefactors (`q_unstated_scaling_prefactors_and_alpha_fit`):*
    - Curve reproductions at printed $\alpha = 0.032$, $\beta = 2.75\times 10^{-4}$, $\epsilon = 10\,k_B T$:
      - `cmp_koff_arrhenius_vs_fig2c_curve`: **reproduced (13/13)**, implied prefactor $k_{\text{off}}(0)/D = 0.1003$ ($10 k_{\text{off}}(0)/D = 1.003$).
      - `cmp_kon_constant_vs_fig2c_curve`: **reproduced (8/8)**, implied prefactor `0.09122` $\times\,5\,\mu\text{M} = 0.4561$.
      - `cmp_kon_adsorption_vs_fig2c_curve`: **close (14/15 reproduced, 1/15 close at $\phi = 0.36$: `0.4728` vs. `[0.4748, 0.4897]`)**, implied prefactor `0.1379` $\times\,5\,\mu\text{M} = 0.6895$.
      - `cmp_kon_coat_hindered_vs_fig2c_curve`: **reproduced (12/12)**, implied prefactor `0.2107` $\times\,5\,\mu\text{M} = 1.054$.
    - Fit quality against the 5 digitized points **with their $\pm 20\%$ error bars**:
      - `cmp_koff_arrhenius_vs_fig2c_data`: **reproduced (5/5)** (all 5 points fall inside their $\pm 20\%$ bars, with $\phi = 6$ just inside at `0.01480` vs. `[0.01443, 0.02311]`).
      - `cmp_kon_coat_hindered_vs_fig2c_data`: **close (4/5 reproduced, 1/5 close at $\phi = 6$: `0.1536` vs. `[0.1688, 0.2627]`)**.
      - `cmp_kon_constant_vs_fig2c_data` and `cmp_kon_adsorption_vs_fig2c_data`: both **not reproduced** (`2 reprod, 2 close, 1 not reprod` at $\phi = 6$).

- **Fig. 2d (p. 4 / 450) — Read (`fig2d_npeak_data`, `fig2d_npeak_constant_curve`, `fig2d_npeak_adsorption_curve`, `fig2d_npeak_coat_hindered_curve`):**
  - *Series:* Peak filament number $N_{\text{peak}}$ vs. $\phi \in [0, 10]$:
    1. `fig2d_npeak_data` (`N_peak_fig2d_data`, `origin: derived`, `error_bars: range`): **5** grey-filled black squares with vertical error bars:
       - $\phi = 0.204$: $7.60$ (Table 1: `7.6`, `0.00%`), bar `[4.23, 10.97]` (`-44.3% / +44.3%`)
       - $\phi = 0.504$: $12.77$ (Table 1: `12.8`, `-0.24%`), bar `[7.28, 18.42]` (`-43.0% / +44.3%`)
       - $\phi = 1.004$: $21.08$ (Table 1: `21.1`, `-0.09%`), bar `[11.85, 30.36]` (`-43.8% / +44.0%`)
       - $\phi = 3.003$: $30.20$ (Table 1: `30.2`, `0.00%`), bar `[16.90, 43.51]` (`-44.1% / +44.1%`)
       - $\phi = 5.989$: $32.29$ (Table 1: `32.3`, `-0.05%`), bar `[18.10, 46.47]` (`-43.9% / +43.9%`)
       This settles `q_npeak_error_propagation_unstated`: all 5 bars are **symmetric $\pm 44\%$** ($1.20^2 - 1 = 0.44$ applied symmetrically around $N_{\text{peak}}$).
    2. `fig2d_npeak_constant_curve` (`N_peak_fig2d_constant`, `origin: predicted`): solid black curve sampled at **10** $\phi$ positions (`0.65..3.45`, $N_{\text{peak}} = 7.92\text{–}47.51$).
    3. `fig2d_npeak_adsorption_curve` (`N_peak_fig2d_adsorption`, `origin: predicted`): dashed black curve sampled on **9** dashes (`0.70..2.40`, $N_{\text{peak}} = 12.16\text{–}48.96$).
    4. `fig2d_npeak_coat_hindered_curve` (`N_peak_fig2d_coat_hindered`, `origin: predicted`): dotted black saturating curve (misreferenced as `Fig. 3d` on p. 453) sampled at **13** $\phi$ positions (`0.35..9.80`, $N_{\text{peak}} = 13.13\text{–}26.92$).
  - *Calibration (1200 dpi):* Linear $x$-axis from left/right spines + 9 inward ticks (`0..10`, 11 marks, residual `3.54 px` across `2544 px` = `0.14%`); linear $y$-axis from bottom spine + 9 clean inward ticks (`0, 10, 15, ..., 50`, 10 marks, residual `1.94 px` across `2496 px` = `0.08%`; tick `5` occluded by starting curves).
  - *Verdicts & cross-panel consistency:*
    - Curve reproductions at printed $\alpha = 0.032$, $\beta = 2.75\times 10^{-4}$, $\epsilon = 10\,k_B T$:
      - `cmp_N_peak_constant_vs_fig2d_curve`: **reproduced (10/10)**, implied prefactor `0.2096` $\times\,25\,\mu\text{M}^2 = 5.240$ (matches Fig. 2c's $(0.4561 / (2 \times 0.1003))^2 = 5.17$).
      - `cmp_N_peak_adsorption_vs_fig2d_curve`: **reproduced (9/9)**, implied prefactor `0.4755` $\times\,25\,\mu\text{M}^2 = 11.89$ (matches Fig. 2c's $(0.6895 / (2 \times 0.1003))^2 = 11.81$).
      - `cmp_N_peak_coat_hindered_vs_fig2d_curve`: **close (12/13 reproduced, 1/13 close at $\phi = 0.85$: `19.58` vs. `[19.66, 20.22]`)**, implied prefactor `1.117` $\times\,25\,\mu\text{M}^2 = 27.93$ (matches Fig. 2c's $(1.054 / (2 \times 0.1003))^2 = 27.61$).
    - Fit quality against the 5 digitized points **with their symmetric $\pm 44\%$ error bars**:
      - `cmp_N_peak_coat_hindered_vs_fig2d_data`: **reproduced (5/5)** — every single fitted $N_{\text{peak}}$ falls cleanly inside its $\pm 44\%$ error bar.
      - `cmp_N_peak_constant_vs_fig2d_data` (`2 reprod, 2 close, 1 not reprod`) and `cmp_N_peak_adsorption_vs_fig2d_data` (`4 close, 1 not reprod`): both **not reproduced**, diverging exponentially past $\phi \ge 3$.

- **Fig. 3a–b (p. 5 / 451) — Left out:**
  - *Reason:* Schematic cartoons of disc-like bundle nucleus growth dynamics showing short actin cylinders and oval fascin cross-linkers at low (a) and high (b, fascin-coat hindrance) $[F]/[A]$ ratios; drawings with no numeric axes.

- **Fig. 4a–d (p. 6 / 452) — Left out:**
  - *Reason:* Fluorescence microscopy images (a, b; scale bars $20\,\mu\text{m}$) and cryo-TEM micrographs (c, d; scale bars $0.2\,\mu\text{m}$) comparing preformed F-actin + fascin (a, c) vs. monomeric co-assembly (b, d) at $[A]/[F] = 1/1$, $[A] = 5\,\mu\text{M}$; micrographs with no numeric axes.

- **Pages 7–8 (453–454) — No figures:** Continuation of results, discussion, acknowledgments, and references (text only).

- **Fig. S1a–h (caption on p. 9 / SI-1) — Left out (images missing from PDF):**
  - *Reason:* Supporting Information p. 9 prints the caption for Figure S1 describing linear (`a, c, e, g`) and semi-log (`b, d, f, h`) plots of $P(N_b)$ and Eq. (1) fits for $[F]/[A] = 1/5, 1/2, 1/1, 6/1$, but the actual plot image panels are absent from `2008_Haviv.pdf` (pages 9–10 converted from Springer file `249_2007_236_MOESM1_ESM.doc` contain only captions on p. 9 and Table S1 on p. 10; logged in `q_missing_fig_s1_plots_in_pdf`).

- **Page 10 (SI-2) — No figures:** Supporting Table S1 only (all tabular values already recorded in `data.yaml`).
