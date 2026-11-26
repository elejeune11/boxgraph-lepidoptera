# Figures in Henle & Pincus (2005) (`henle2005`)

All four pages of `corpus/actin_bundling/2005_Henle.pdf` (`pages/page-001.png`–`page-004.png`, journal pages `060801-1` to `060801-4`) were inspected. Pages 1, 2, and 4 contain only text, equations, and references (no figures or tables). Page 3 contains two theoretical figures, each with a main panel and an inset (four plotted panels total, 11 distinct curves/lines, 196 digitized points). The paper contains no experimental measurements, error bars, bar charts, or histograms. Every figure panel and every plotted curve/line is digitized with calibrated `mode: code` readers in `figures/` (shared helper `figures/_henle_figs.py`, using the passed image array directly via `scipy.ndimage.find_objects`, ~230–330 ms per reader) and saved in `readings.json`:

- **Fig. 1 main (p. 3, 1200 dpi, crop `[1200, 880, 4560, 2680]`) — read (`fig1_nq_01uM`, `fig1_nq_10uM`, `fig1_nq_1mM`; 61 points total)**:
  - Curves: $\beta\Delta F/(NL)$ vs $R/a$ for pointlike counterions at $n_q = 0.1\,\mu\text{M}$ (`fig1_nq_01uM`, solid curve, 10 sampled points at $R/a \in [2.5, 9.5]$), $10\,\mu\text{M}$ (`fig1_nq_10uM`, dashed curve, 25 dash midpoints over $R/a \in [1.23, 9.70]$), and $1\,\text{mM}$ (`fig1_nq_1mM`, dotted curve, 26 dot centroids over $R/a \in [1.30, 9.69]$).
  - Calibration & misfit: frame borders and inward ticks matched to 4 major x-ticks ($R/a = 2, 4, 6, 8$; $278.3\,\text{px/unit}$, max linear-fit misfit $0.80\,\text{px} = 0.0029$) and 13 regular y-ticks ($2.25$ down to $-0.75\,\text{nm}^{-1}$ every $0.25\,\text{nm}^{-1}$; $440.8\,\text{px/(nm}^{-1})$, max misfit $2.10\,\text{px} = 0.0048\,\text{nm}^{-1}$), giving calibrated frame $[0.999, 9.982] \times [-1.002, 2.501]$.
  - Note: the rectangular inset exclusion mask skips three dots of `fig1_nq_1mM` near $R/a \approx 4.90, 5.20, 5.50$ (preserved unchanged from the audited first-reader digitization; all 26 retained dots sit dead-center on their printed markers).

- **Fig. 1 inset (p. 3, 2400 dpi, crop `[5600, 1880, 8960, 3680]`) — read (`fig1_inset_lam_tot`, `fig1_inset_lam_tot_asymp`; 22 points total)**:
  - Curves: total effective linear charge density $\lambda_{\text{tot}}(R)$ at $n_q = 10\,\mu\text{M}$ (`fig1_inset_lam_tot`, solid curve, 12 sampled points at $R/a \in [5, 95]$) and its horizontal large-$R$ asymptotic line (`fig1_inset_lam_tot_asymp`, dashed line, 10 interior dash centroids over $R/a \in [11.83, 90.25]$, reading $8.195 \pm 0.088\,\text{nm}^{-1}$).
  - Calibration & misfit: frame borders and inward ticks matched to 9 regular bottom-spine x-ticks ($R/a = 10, 20, \dots, 90$; $22.6\,\text{px/unit}$, max misfit $3.90\,\text{px} = 0.173$) and 6 regular left-spine y-ticks ($\lambda_{\text{tot}} = 8, 7, 6, 5, 4, 3\,\text{nm}^{-1}$; $197.6\,\text{px/(nm}^{-1})$, max misfit $3.75\,\text{px} = 0.019\,\text{nm}^{-1}$), giving calibrated frame $[0.970, 100.171] \times [1.576, 8.580]$.

- **Fig. 2 main (p. 3, 1200 dpi, crop `[5600, 880, 8960, 2680]`) — read (`fig2_nq_01uM`, `fig2_nq_10uM`, `fig2_nq_1mM`; 67 points total)**:
  - Curves: $\beta\Delta F/(NL)$ vs $R/a$ for finite-size counterions ($\gamma = 1.7\,\text{nm}^{-1}$, $\lambda^* = 0.01\lambda_0$) at $n_q = 0.1\,\mu\text{M}$ (`fig2_nq_01uM`, solid curve, 10 sampled points at $R/a \in [7, 28]$), $10\,\mu\text{M}$ (`fig2_nq_10uM`, dashed curve, 29 dash midpoints over $R/a \in [1.53, 28.69]$), and $1\,\text{mM}$ (`fig2_nq_1mM`, dotted curve, 28 dot centroids over $R/a \in [3.73, 27.99]$, bounded at $x_{\max} = 29.2$ to avoid border crowding).
  - Calibration & misfit: frame borders and inward ticks matched to 5 regular bottom-spine x-ticks ($R/a = 5, 10, 15, 20, 25$; $86.0\,\text{px/unit}$, max misfit $1.20\,\text{px} = 0.014$) and 9 regular left-spine y-ticks ($-5.1$ down to $-5.9\,\text{nm}^{-1}$ every $0.1\,\text{nm}^{-1}$; $1539.0\,\text{px/(nm}^{-1})$, max misfit $1.41\,\text{px} = 0.00092\,\text{nm}^{-1}$), giving calibrated frame $[1.009, 29.986] \times [-5.998, -5.001]$.

- **Fig. 2 inset (p. 3, 2400 dpi, crop `[14320, 1880, 17760, 3680]`) — read (`fig2_inset_nq_01uM`, `fig2_inset_nq_10uM`, `fig2_inset_nq_1mM`; 46 points total)**:
  - Curves: $\beta\Delta F/(NL)$ vs $R/a$ for frustrated attractive interactions ($\xi = 0.1$, $\phi_{\max} = 2$, $\gamma = 1.7\,\text{nm}^{-1}$) at $n_q = 0.1\,\mu\text{M}$ (`fig2_inset_nq_01uM`, solid curve, 9 sampled points at $R/a \in [5.5, 13]$), $10\,\mu\text{M}$ (`fig2_inset_nq_10uM`, dashed curve, 7 interior dash midpoints over $R/a \in [4.42, 13.31]$, excluding the bottom-clipped dash near $R/a \approx 3.6$), and $1\,\text{mM}$ (`fig2_inset_nq_1mM`, dotted curve, 30 dot centroids over $R/a \in [1.20, 13.80]$).
  - Calibration & misfit: frame borders and inward ticks matched to 6 regular top-spine x-ticks ($R/a = 2, 4, 6, 8, 10, 12$; $172.5\,\text{px/unit}$, max misfit $3.54\,\text{px} = 0.021$) and 10 regular left-spine y-ticks ($0.3$ down to $-0.6\,\text{nm}^{-1}$ every $0.1\,\text{nm}^{-1}$; $1258.2\,\text{px/(nm}^{-1})$, max misfit $3.42\,\text{px} = 0.0027\,\text{nm}^{-1}$), giving calibrated frame $[1.016, 14.004] \times [-0.699, 0.399]$.

## Left Out

No figures, panels, or curves were left out (all 4 panels and all 11 plotted curves/lines on p. 3 are digitized; pp. 1, 2, and 4 have no figures).

## Relevance to Open Questions

- `q_n1_10nM_vs_10mM`: All 8 digitized curves across Fig. 1 main, Fig. 1 inset, and Fig. 2 main fail reproduction at the caption's literal `n_1 = 10 nM`, and all 8 simultaneously reproduce within reading uncertainty under the diagnostic `n_1 = 10 mM` unit-typo hypothesis (including the Fig. 1 inset asymptotic line: computed $8.190\,\text{nm}^{-1}$ vs digitized $8.195 \pm 0.088\,\text{nm}^{-1}$).
- `q_finite_size_eq6_rhs_sign`: The three digitized curves of Fig. 2 main confirm that direct minimization of the printed finite-size free energy (`finite_size_free_energy_p3`) and the corrected LHS Eq. (6) variant (`finite_size_eq6_corrected`) reproduce the plotted curves under `n_1 = 10 mM`, whereas the literal printed RHS sign has no real solution.
- `q_fig2_inset_unexplained_discrepancy`: Even under `n_1 = 10 mM`, the 46 digitized points of Fig. 2 inset are not reproduced by the printed frustrated-attraction model (`frustrated_bundle_p4`), confirming a persistent quantitative discrepancy beyond the `n_1` typo.

## Overlay Verification

All 11 overlays were generated via `boxgraph figure check` into `figwork/<reading_id>.png`, assembled into `figwork/contact_sheet_fig1.png` and `figwork/contact_sheet_fig2.png`, and visually inspected: all 196 magenta crosses sit cleanly on their target curves, dash midpoints, and dot centroids with monotonic `x`, zero missing values, and verified linear calibration misfits.
