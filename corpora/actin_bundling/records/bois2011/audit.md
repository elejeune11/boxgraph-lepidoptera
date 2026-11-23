# Audit Log — `bois2011`

Independent audit of `corpora/actin_bundling/records/bois2011/` against `pages/page-001.png`–`pages/page-004.png` and `text/page-001.txt`–`text/page-004.txt` (Bois, Jülicher & Grill, *Phys. Rev. Lett.* **106**, 028103, 2011).

## Summary of verification

- **Data transcription check (`data.yaml`, `source.yaml`, `references.yaml`)**: All 22 original data entries, bibliographic metadata, and 12 references were checked against the page images (`page-001.png`–`page-004.png`). No digit, sign, exponent, unit, page, or quote transcription errors were found (0 transcription errors fixed).
- **Model equation check (`models/*`)**: All 23 original models were checked against printed Eqs. (1)–(15) and inline formulas on pp. 1–4, and all 18 `symbolic()` definitions were verified against `run()`. No equation, sign, or factor transcription errors were found (0 equation/sign/factor errors fixed).

## Changes made (missed printed data, model, run, and comparison)

- Added quantities `regulator_concentration` and `concentration_gradient` in `quantities.yaml` (p. 1, `028103-1`, Eqs. 1 and 5) to support the omitted no-flux boundary condition derivation from Eq. (5).
- Added missed printed boundary data `v_noflux_boundary` (`0 m/s`), `dc_dx_noflux_boundary` (`0 1/m`), and missed printed derived result `j0_steady_noflux` (`0 m/s`) in `data.yaml` (p. 1, `028103-1`, text after Eq. 5: *"For no-flux boundary conditions, $v(0) = \partial_x c(0) = 0$ requires that $j_0 = 0$"*).
- Added unstated boundary concentration `c_boundary_noflux_added` (`1`, alternatives `["0.5", "2"]`) in `added.yaml` (p. 1, `028103-1`, Eq. 5).
- Added missed model `steady_flux_eq5` (`models/steady_flux_eq5/model.yaml` and `model.py`) implementing Eq. (5) (`v = D * d_x(ln c) + j0 / c`, solved for `j0 = c * v - D * d_x c`) (p. 1, `028103-1`, Eq. 5).
- Added run `run_steady_flux_noflux` in `runs.yaml` and reproduction comparison `cmp_j0_steady_noflux` in `comparisons.yaml` checking `j0_noflux_pred` against `j0_steady_noflux` (p. 1, `028103-1`, Eq. 5).

## Disagreements with first-reader choices (logged as `reading` questions)

- `q_audit_asdm_pde_grid_negativity` (p. 4, `028103-4`, Fig. 5(d)): Disagreement with choosing `n_grid: 64` as the default grid in `asdm_pde_eq14` and refining only `a_max` in `ref_asdm_pde_fig5d_ngrid`; at `n_grid = 64`, `run_asdm_pde_fig5d` outputs an unphysically negative concentration `a_min_fig5d_pred = -0.00575` (due to spectral Gibbs undershoot, and 2/3 dealiasing does not fully dealias the cubic reaction term `a^2 * s`) and an 8.2% error in `s_min` (`0.09142` vs. `0.09959`/`0.09956` at `n_grid = 128`/`256`, where `a_min` converges to `+0.00273`/`+0.00268`).
- `q_audit_pde_hardcoded_initial_conditions` (pp. 3–4, `028103-3`–`028103-4`, Figs. 4 and 5(d)): Disagreement with hardcoding critical initial-perturbation modal shapes inside `single_species_pde_eq1_3/model.py` (mode-1 symmetry-breaking ratio `0.1`, where lowering to `1e-5` flips `cmp_final_peak_count_fig4` from `1` to `2` peaks at `T_end = 20`) and `asdm_pde_eq14/model.py` (pure 3-mode perturbation enforcing 3-fold symmetry at `T_end = 40`, whereas mode `n = 4` has the largest linear growth rate by Eq. 15 and remains 4-peaked at `T_end = 40` before coarsening to 3 peaks by `T ~ 200`) instead of exposing them in `added.yaml` with alternatives.
- `q_audit_exact_bc_added_alternatives` (pp. 1–2, `028103-1`–`028103-2`, Eqs. 7–8): Disagreement with entering exact mathematical consequences of the paper's stated boundary and homogeneity conditions (`delta_sigma_periodic_added`, `c_ratio_periodic_added`, `du_dx_homog_added`) in `added.yaml` with nonzero/non-unit alternatives that violate periodicity or homogeneity and falsely make exact identities (`cmp_j0_steady_periodic`, `cmp_v_homog_steady`) look sensitive to unstated inputs.

## New open questions added (non-`reading`)

- `q_audit_fig5a_L_over_l_irrelevant` (`possible_error`, p. 4, `028103-4`, Fig. 5(a)): Fig. 5(a)'s caption specifies "assuming $L/\ell \gg 1$" for the `Pe = 0` stability diagram, but when `Pe = 0` the hydrodynamic length $\ell = \sqrt{\eta/\gamma}$ does not appear in Eq. (15); the relevant continuum-wavenumber condition for Fig. 5(a) is $L \gg \sqrt{D_s/\rho_s}$, whereas $L/\ell \gg 1$ applies to Fig. 5(b) ($\mathrm{Pe} > 0$).
