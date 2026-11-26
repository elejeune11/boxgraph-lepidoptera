# Figure Inspection — `gov2008`

Direct visual inspection of all five pages of Nir S. Gov, *Phys. Rev. E* **78**, 011916 (2008), "Packing defects and the width of biopolymer bundles" (`corpus/actin_bundling/2008_Gov.pdf`, PDF pages 1–5, journal pages `011916-1`–`011916-5`). All five page PNGs (`pages/page-001.png`–`pages/page-005.png`) were inspected directly, and the two figure-bearing pages (pages 1 and 3) were additionally rendered at 300 dpi with pixel rulers (`figwork/p1.png`, `figwork/p3.png`) via `boxgraph figure page` and cropped at **900 dpi** (`figwork/fig1_900.png`, `figwork/fig2_900.png`) for close inspection of every symbol, arrow, and label.

**Summary:** This is a purely analytical paper containing only two figures (**Fig. 1** on page 1 and **Fig. 2** on page 3, 6 panels total) and **no tables**. Every figure panel is a qualitative geometric schematic with symbolic labels only — no numeric axes, tick marks, scale bars, plotted model curves, measurements, error bars, or histograms appear anywhere in the paper. Checking `data.yaml` also confirms that zero existing data entries come from figures (the only approximate `"~"` entry, `filament_radius_a_printed = "~10"` nm, is cited from prose on p. 2 col 1). Accordingly, **0 calibrated figure readings** were created and `figures.yaml` is omitted.

---

## Page-by-page and panel-by-panel inventory

### Page 1 (`011916-1`, `pages/page-001.png`, `figwork/p1.png`, `figwork/fig1_900.png`) — Fig. 1 (3 panels)

- **Fig. 1(a) — Left out (purely schematic drawing):**
  - *Visual contents:* 3D schematic cylinder (dashed outline) of a biopolymer bundle of radius $R$ (horizontal single-headed arrow from cylinder center to right rim) and overall length $L_b$ (vertical double-headed arrow along the right side), containing straight white cylindrical filaments and colored (red/orange) defective filaments. Depicts a short interrupted filament near the top left (splay defect) and pairs of crossing/entangled filaments (twist defects) bending over two horizontal dotted defect planes, each marked by a vertical double-headed arrow labeled $L_p$ (healing length).
  - *Reason left out:* Qualitative geometric schematic with symbolic labels ($R$, $L_b$, $L_p$) only; no numeric axes, ticks, scale bars, or plotted numerical data. (The caption states the symbolic curvature scaling $\sim a/L_p^2$ [14], which is already captured analytically in the record's models rather than plotted.)

- **Fig. 1(b) — Left out (purely schematic drawing; corroborates `q_audit_filament_radius_vs_spacing_ambiguity`):**
  - *Visual contents:* 2D cross-sectional schematic far below or above a defect plane showing 19 circular filaments arranged on a regular triangular/hexagonal lattice inside a dashed bundle outline, with two central filaments shaded orange-red, linked by a horizontal dotted segment, and accompanied by a curved rotation arrow indicating a twist defect pair. A single horizontal double-headed arrow labeled $a$ spans the **center-to-center distance** between two adjacent filaments in the top row of the lattice.
  - *Reason left out:* Qualitative cross-sectional schematic with a single symbolic spacing label $a$; no numeric axes, ticks, scale bar, or numerical values.
  - *Relation to open questions:* High-resolution (900 dpi) inspection directly confirms the observation in `q_audit_filament_radius_vs_spacing_ambiguity`: Fig. 1(b) unmistakably draws $a$ as the center-to-center lattice spacing (effective filament diameter), whereas p. 2 col 1 prose describes $a \sim 10\text{ nm}$ as "the radius of the individual filaments and surrounding cross linkers [5,6]" (matching the area ratio $(R/a)^2$ in Eqs. 1–6). Because Fig. 1(b) has no numeric scale bar, it confirms the schematic-vs.-prose ambiguity without numerically settling whether $10\text{ nm}$ was intended as radius or center-to-center spacing.

- **Fig. 1(c) — Left out (purely schematic drawing; relates to `q_continuum_elasticity_multipole_decay_2d_vs_3d`):**
  - *Visual contents:* 2D cross-sectional schematic in the defect plane showing the orange-red twist pair rotated vertically (connected by a vertical dotted line), displacing the surrounding 17 filaments (solid circles) away from their undistorted hexagonal lattice sites (dashed circles). Qualitative peripheral arrows indicate inward compression along the horizontal axis (2 inward arrows, left and right) and outward expansion along the vertical and diagonal directions (6 outward arrows).
  - *Reason left out:* Qualitative displacement/strain-field schematic with no axes, ticks, scale bar, or quantitative displacement magnitudes.
  - *Relation to open questions:* Visually matches the quadrupolar symmetry described on p. 2 col 2, but carries no quantitative radial profile and therefore cannot settle the 2D ($1/r^3$) vs. 3D ($1/r^4$) quadrupole decay discrepancy noted in `q_continuum_elasticity_multipole_decay_2d_vs_3d`.

### Page 2 (`011916-2`, `pages/page-002.png`) — No figures

- **Page 2 — No figures or tables:** Two columns of prose and Eq. (1). No figure panels or tables appear on this page.

### Page 3 (`011916-3`, `pages/page-003.png`, `figwork/p3.png`, `figwork/fig2_900.png`) — Fig. 2 (3 panels)

- **Fig. 2(a) — Left out (purely schematic drawing):**
  - *Visual contents:* 3D schematic cylinder of radius $R$ (horizontal arrow from center to right rim) illustrating the empty-shell model (Eqs. 2–3): an outer concentric ring of thickness $\lambda$ (double-headed arrow across the outer rim) is defect-free, while the inner circular core contains 4 red dumbbell pairs representing localized packing defects with typical inter-defect spacing $L$ (double-headed arrow between two pairs, with leader line to label $L$; caption gives areal density $\rho = 1/L^2$ and core energy $E_c$).
  - *Reason left out:* Qualitative model schematic with symbolic labels ($R$, $\lambda$, $L$) only; no numeric axes, ticks, scale bar, or plotted data.

- **Fig. 2(b) — Left out (purely schematic drawing; relates to `q_eqs7_to_9_assumption_Ec_lt_Ecstar_vs_lambda_zero`):**
  - *Visual contents:* 3D schematic cylinder of radius $R$ illustrating uniform bulk monopole defects with long-range logarithmic elastic repulsion (Eq. 4): 9 red dumbbell pairs distributed uniformly across the entire top circular cross-section (with spacing arrow labeled $L$, no empty outer shell) and 7 `+` symbols on the bottom cross-section depicting the analogous 2D uniform gas of same-sign logarithmic charges.
  - *Reason left out:* Qualitative conceptual schematic with symbolic labels ($R$, $L$, `+`) only; no numeric axes, ticks, scale bar, or plotted data.
  - *Relation to open questions:* Consistent with `q_eqs7_to_9_assumption_Ec_lt_Ecstar_vs_lambda_zero`, Fig. 2(b) draws defects extending all the way across the full cross-section with no defect-free outer shell ($\lambda = 0$), matching the absence of Eq. (2)'s shell term in Eqs. (7) and (9).

- **Fig. 2(c) — Left out (purely schematic drawing; relates to `q_eqs7_to_9_assumption_Ec_lt_Ecstar_vs_lambda_zero`):**
  - *Visual contents:* 3D schematic cylinder of radius $R$ illustrating surface monopole defects (Eq. 5): 6 red dumbbell pairs distributed along the perimeter of the top cross-section (with spacing arrow labeled $L$) and 6 `+` symbols along the perimeter of the bottom cross-section depicting same-sign logarithmic charges confined to the bundle surface.
  - *Reason left out:* Qualitative conceptual schematic with symbolic labels ($R$, $L$, `+`) only; no numeric axes, ticks, scale bar, or plotted data.

### Page 4 (`011916-4`, `pages/page-004.png`) — No figures

- **Page 4 — No figures or tables:** Two columns of prose and Eqs. (9)–(11). No figure panels or tables appear on this page.

### Page 5 (`011916-5`, `pages/page-005.png`) — No figures

- **Page 5 — No figures or tables:** Conclusions, Acknowledgments, and References [1]–[24]. No figure panels or tables appear on this page.
