# Cross-Record `relate` & `validate` Audit — `kumar2014` (Pulsatory Patterns in Active Fluids)

Owner record: `kumar2014` (K. Vijay Kumar, J. S. Bois, F. Jülicher, and S. W. Grill, *Phys. Rev. Lett.* **112**, 208101, 2014).
Audited ready pairs (all 5 active-gel mechanochemical patterning pairs):
1. `bois2011 + kumar2014` (direct precursor; explicitly cited as `ref-10` on pp. 1–3)
2. `hannezo2015 + kumar2014` (independent active-gel-with-turnover extension of `bois2011`)
3. `kumar2014 + mietke2019` (deformable zero-friction active surfaces; `mietke2019` cites `kumar2014` as Ref. 15 on pp. 1, 3, 5)
4. `barberi2023 + kumar2014` (localized snaking active-gel contractions; `barberi2023` cites `kumar2014` as Ref. [13] on p. 1)
5. `kumar2014 + mirza2025` (compressible active-nematic gels; `mirza2025` cites `kumar2014` on pp. 2, 3, 8)

---

## Summary Counts (inside `kumar2014/`)

- **`cited.yaml`**: **0** links (`kumar2014` is a purely dimensionless theoretical study; its sole `from_source` datum in `data.yaml` is `eta_v_over_eta_2d = 3` from `ref-29`, Jenkins & Dysthe 1997, outside the corpus; zero empirical or numerical values are imported from `bois2011`).
  - *Read-only partner reference links noted (not edited here per strict ownership rule)*: `mietke2019/references.yaml` (Ref. 15), `barberi2023/references.yaml` (Ref. [13]), and `mirza2025/references.yaml` (`ref-kumar2014`) all cite `kumar2014` in their bibliographies but currently lack `source: kumar2014` because `kumar2014` was ingested after those records.
- **`relations.yaml`**: **42** cross-record relations added (**15** `bois2011`, **8** `hannezo2015`, **6** `mietke2019`, **5** `barberi2023`, **8** `mirza2025`), alongside the 6 existing within-record relations (**48/48 `verified`**, 0 `ambiguous`, 0 `not_verified`, 0 `refuted`, 0 `unverifiable`).
- **`validations.yaml`**: **0** entries (empty list `[]` with explicit audit header; 0 genuine non-forced cross-validations across all 5 pairs).
- **`questions.yaml`**: **4** cross-paper comparison & audit questions appended (`q_cross_bois2011_single_regulator_and_asdm_contrast`, `q_cross_viscous_conventions_and_2d_factor4`, `q_cross_ard_hopf_vs_single_species_turnover_and_oscillation_mechanisms`, `q_cross_validation_audit_all_5_partners`).

---

## Pair-by-Pair Audit

### 1. `bois2011 + kumar2014` (Direct 2-Species Extension of Bois et al. 2011)
- **Citations / Relations / Validations**: `ref-10` (`source: bois2011`) in `references.yaml`, 0 `cited.yaml` numeric links / **15** verified relations (**8 `stated`**, **7 `proposed`**) / 0 validations (both are pure theory/simulation with 0 `origin: measured` datums).
- **Verified relations (15/15)**:
  - *Shared active-gel definitions (`same_as`, `stated` p. 2)*:
    - `rel_bois2011_hydrodynamic_length`: `hydrodynamic_length` $\ell = \sqrt{\eta/\gamma}$ (symbolic).
    - `rel_bois2011_characteristic_velocity`: `characteristic_velocity` $U = (\zeta\Delta\mu)_0/\sqrt{\eta\gamma}$ (symbolic).
    - `rel_bois2011_diffusive_timescale`: `diffusive_timescale` $\tau = \ell^2/D$ (symbolic).
    - `rel_bois2011_peclet_number_physical`: `peclet_number_physical` $P = (\zeta\Delta\mu)_0/(\gamma D)$ (symbolic).
    - `rel_bois2011_peclet_number_timescale`: `bois2011/peclet_number_timescale` ($U\ell/D$) ↔ `kumar2014/peclet_number_timescale` ($\tau/\tau_a$) (symbolic).
    - `rel_bois2011_discrete_wavenumber`: `discrete_wavenumber` $k_n\ell = \omega n\pi/(L/\ell)$ (symbolic; Bois's BC factor $\beta \in \{1,2\}$ is renamed $\omega$ in Kumar because $\beta$ is reassigned to active-stress asymmetry in Eq. 11).
  - *Exact 1-species reductions of linear stability & thresholds*:
    - `rel_bois2011_eq7_equal_diff_reduces_eq12` (`reduces_to`, `stated` p. 3, err `8.3e-09`): `linear_stability_ad_eq7` → `bois2011/single_species_dispersion_eq12` as $\alpha \to 1, \nu \to 1$ with effective sensitivity $c_0 f'(c_0) = A_0 f_A + I_0 f_I$.
    - `rel_bois2011_eq7_single_regulator_reduces_eq12` (`reduces_to`, `proposed`, err `6.8e-09`): `linear_stability_ad_eq7` → `bois2011/single_species_dispersion_eq12` as $I_0 f_I \to 0, \nu \to 1$ ($\alpha > 1$).
    - `rel_bois2011_fig1_beta1_reduces_hill` (`reduces_to`, `proposed`, err `1.8e-07`): `linear_stability_ad_fig1` → `bois2011/single_species_dispersion_hill` at $\beta \to 1, \nu \to 1$ with $c_0 = 1/\xi = 1/3$ and $\mathrm{Pe}_{\text{Bois}} = 2 P (1+\xi)^2/\xi = (32/3)P$.
    - `rel_bois2011_fig1d_stationary_same_as_eq13` (`same_as`, `stated` p. 3, max rel diff `2.84e-16`): `linear_stability_boundary_ad_fig1d` (`P_c`) ↔ `bois2011/single_species_instability_eq13` (`Pe_crit`) across the entire 1D ($\nu=1$) two-activator regime $\beta \in (-1, 1)$ with diffusivity-weighted sensitivity $c_0 f'(c_0) = (1+\beta) + (1-\beta)/\alpha$.
    - `rel_bois2011_fig1d_beta1_reduces_eq13` (`reduces_to`, `proposed`, err `6.7e-09`): `linear_stability_boundary_ad_fig1d` → `bois2011/single_species_instability_eq13` at $\beta \to 1, \nu \to 1$ ($\alpha$ cancels identically, leaving $c_0 f'(c_0) = 2$).
    - `rel_bois2011_eq10_reduces_eq13` (`reduces_to`, `proposed`, err `1.5e-08`): `oscillatory_threshold_eq10` → `bois2011/single_species_instability_eq13` as $\alpha \to 0, I_0 f_I \to 0, \nu \to 1$.
    - `rel_bois2011_trace_eq8_reduces_eq12_dim` (`reduces_to`, `proposed`, err `3.2e-07`): dimensional `trace_discriminant_eq8_9` (`tr_L` [1/s]) → `bois2011/single_species_dispersion_eq12` (`lam` [1/s]) as $\alpha \to 0, I_0 f_I \to 0, \nu \to 1$.
  - *Two-species reaction-diffusion-advection bridges*:
    - `rel_bois2011_asdm_zero_rxn_reduces_eq7` (`reduces_to`, `proposed`, err `1.5e-07`): `bois2011/asdm_dispersion_eq15` → `kumar2014/linear_stability_ad_eq7` as $\rho_a, \rho_s \to 0$ (passive substrate $s \leftrightarrow A$ with $A_0 f_A = 0$; stress-regulating activator $a \leftrightarrow I$ with $I_0 f_I = 1/4, \alpha = D_a/D_s, \nu = 1$).
    - `rel_bois2011_ard_eq12_13_scales_rda_eq15` (`scales_as`, `proposed`, max slope diff `2.4e-11`): `ard_linear_stability_eq12_13` → `bois2011/single_species_rda_dispersion_eq15` at $\rho = 1, \nu = 1$ with $\partial_c R = -\kappa\tau$.

### 2. `hannezo2015 + kumar2014` (Active-Gel Turnover & Patterning)
- **Citations / Relations / Validations**: 0 (independent extensions of `bois2011`) / **8** verified relations / **0** validations (`hannezo2015`'s 35 measured datums are stationary single-species *Drosophila* tracheal actin rings, FRAP $D, \tau$, fluorescence ramps, and tube diameters; `kumar2014` has no dimensional stationary wavelength model and no identified second inhibitor species $I$).
- **Verified relations (8/8)**:
  - *Shared scales (`same_as`, symbolic)*: `rel_hannezo2015_hydrodynamic_length` ($\ell = l_h$), `rel_hannezo2015_turnover_length` ($\ell_\kappa = \sqrt{D/\kappa} = \sqrt{D\tau} = l_t$), `rel_hannezo2015_diffusive_timescale` ($T = \tau$), `rel_hannezo2015_peclet_contractility` ($P = \chi$ under thin-film thickness cancellation).
  - *Exact 1-species ARD dispersion reductions (`reduces_to`)*: `rel_hannezo2015_ard_equal_turnover_reduces_eq8` ($\rho \to 1, \nu \to 1$, err `8.1e-09`) and `rel_hannezo2015_ard_single_reg_reduces_eq8` ($I_0 f_I \to 0, \nu \to 1$, $\rho > 1$, err `3.1e-08`) both recover `hannezo2015/dispersion_relation_dimless_eq8` identically.
  - *Exact Hopf-to-stationary threshold identity*:
    - `rel_hannezo2015_ard_osc_same_as_eq10` (`same_as`, max rel diff `0.00e+00`): `ard_stability_boundary_2d` (`P_c_osc`) ↔ `hannezo2015/dimensionless_critical_contractility_eq10` (`chi_c`) with $\phi = \frac{1}{2}\nu\kappa\tau(1+\rho)$ across all $(\beta, \kappa\tau, \rho, \nu) > 0$, because the factor of 2 in $\operatorname{tr}(\tau L_{\text{ARD}}) = 0$ cancels $A_0 f_A + I_0 f_I = 2$.
    - `rel_hannezo2015_ard_stat_scales_eq10` (`scales_as`, max slope diff `4.4e-13`): `ard_stability_boundary_2d` (`P_c_stat`) → `hannezo2015/dimensionless_critical_contractility_eq10` at $\beta = 1$ ($\phi = \nu\kappa\tau$, exact constant factor $1/2$ from $A_0 f_A = 2$).

### 3. `kumar2014 + mietke2019` (Fixed-Substrate Two-Species Oscillations vs Deformable Active Surfaces)
- **Citations / Relations / Validations**: `mietke2019` cites `kumar2014` as Ref. 15 (pp. 1, 3, 5) / **6** verified relations / 0 validations (both pure theory/simulation).
- **Verified relations (6/6)**:
  - `rel_mietke2019_eq11_scales_hill_eq12` (`scales_as`, max slope diff `0.0e+00`): `active_stress_regulation_eq11` (at $f_0 = 0, \beta = 1$) → `mietke2019/hill_regulation_eq12` ($m = 1, c_0 = \xi$; exact factor of 2).
  - `rel_mietke2019_nu_bar_scales_viscous_coeff_2d` (`scales_as`, max slope diff `2.2e-13`): `mietke2019/viscosity_parameter_nu_bar_corrected` ($\bar{\nu} = (\eta_s+\eta_b)/(4\eta_s)$) → `kumar2014/viscous_coefficient_2d` ($\nu = 1 + \eta_v/\eta$; exact factor $1/4$).
  - `rel_mietke2019_eq10_reduces_cyl_eq27` (`reduces_to`, err `4.1e-16`): `oscillatory_threshold_eq10` → `mietke2019/cylinder_critical_contractility_eq27` ($k = 0$) as $L/\ell \to 0$ (zero-friction limit).
  - `rel_mietke2019_ard_osc_reduces_cyl_dimless` (`reduces_to`, err `2.2e-07`): `ard_stability_boundary_2d` (`P_c_osc`) → `mietke2019/cylinder_critical_contractility_dimless` ($D_{\text{tilde}} = 0$) as $\kappa\tau \to \infty$ (zero-friction flat limit).
  - `rel_mietke2019_trace_eq8_approx_sphere_polar` (`approximated_by`, max err `7.70e-04 < 0.01`): `trace_discriminant_eq8_9` (`tr_L`) ≈ `mietke2019/sphere_polar_growth_rate_eqs72` (`lambda_1`) in the low-friction ($\nu k^2\ell^2 \gg 1$), single-species ($\alpha = 0, I_0 f_I = 0$) regime.
  - `rel_mietke2019_ard_eq12_13_approx_cyl_jac` (`approximated_by`, max err `8.04e-04 < 0.01`): `ard_linear_stability_eq12_13` (`re_lam_plus`) ≈ `mietke2019/cylinder_jacobian_eqs80_82` (`max_re_lambda_tilde`) at equal turnover ($\rho = 1$), low friction ($\nu k_l^2 \gg 1$), and rigid-cylinder fixed-geometry limit ($\tilde{\kappa} = 10^5$, eliminating shape-driven oscillations).

### 4. `barberi2023 + kumar2014` (Screen-Flow Velocity & 3-Species Excitable vs 2-Species Opposing-Stress Models)
- **Citations / Relations / Validations**: `barberi2023` cites `kumar2014` as Ref. [13] (p. 1) / **5** verified relations / 0 validations (both pure theory/simulation).
- **Verified relations (5/5)**:
  - `rel_barberi2023_hydrodynamic_length` (`same_as`, symbolic): $\ell = \sqrt{\eta/\gamma} = \lambda$.
  - `rel_barberi2023_diffusive_timescale` (`same_as`, symbolic): $\tau = \ell^2/D = \lambda^2/D_i$.
  - `rel_barberi2023_viscous_eq3_derived_nu_2d` (`derived_from`, symbolic `steps: ["subs: d = 2", "simplify"]`): `barberi2023/viscous_stress_1d_eq3_corrected` (`viscous_coeff`) → `kumar2014/viscous_coefficient_2d` ($\nu = 1 + \eta_v/\eta$).
  - `rel_barberi2023_sm15_vel_scales_fourier_eq6` (`scales_as`, max slope diff `2.2e-13`): `barberi2023/linear_stability_1d_sm15` (`vel_response_imag`) → `kumar2014/fourier_velocity_eq6` (`v_hat_amp`) at $B = 0$ (exact log-slope match in both $K$ and $Z$).
  - `rel_barberi2023_sm15_reduces_ard_eq12_13` (`reduces_to`, err `4.9e-09`): `barberi2023/linear_stability_1d_sm15` (`max_growth_rate`) → `kumar2014/ard_linear_stability_eq12_13` (`re_lam_plus`) in the decoupled fast-diffusing nucleator limit ($\Omega_{d,0}, \Omega_d \to 0, D_a \to \infty$, $D_c = 1$).

### 5. `kumar2014 + mirza2025` (Exact 2D $\nu = 4$ Active-Gel & Hexagonal-Droplet Correspondence)
- **Citations / Relations / Validations**: `mirza2025` cites `kumar2014` on pp. 2, 3, and 8 (specifically matching isotropic hexagonal droplets at $\beta = \kappa = \lambda_\odot = 0$ to Kumar 2014 Fig. 3b) / **8** verified relations / 0 validations (both pure theory/simulation).
- **Verified relations (8/8)**:
  - *Shared scales (`same_as`, symbolic)*: `rel_mirza2025_hydrodynamic_length` ($\ell = \ell_s$), `rel_mirza2025_turnover_damkohler_length` ($\ell_\kappa = \ell_D = \sqrt{D/k_d}$), `rel_mirza2025_diffusive_timescale` ($\tau = t_0$).
  - *Exact 2D $\nu = 4$ ($\eta_v = 3\eta$) threshold & wavenumber identities*:
    - `rel_mirza2025_ard_osc_same_as_iso_lam` (`same_as`, max rel diff `0.00e+00`): `ard_stability_boundary_2d` (`P_c_osc`) ↔ `mirza2025/isotropic_threshold_dimless` (`lambda_crit_0`) with $k_d = \nu\kappa\tau(1+\rho)/8$ (at $\nu = 4$, simply the mean turnover rate $k_d = \frac{1+\rho}{2}\kappa\tau$).
    - `rel_mirza2025_ard_osc_k_same_as_iso_nu` (`same_as`, max rel diff `0.00e+00`): `ard_stability_boundary_2d` (`k_c_osc_l`) ↔ `mirza2025/isotropic_threshold_dimless` (`nu_crit_0`) at 2D $\nu = 4$ with $k_d = \frac{1+\rho}{2}\kappa\tau$.
    - `rel_mirza2025_ard_stat_scales_iso_lam` (`scales_as`, max slope diff `4.4e-13`): `ard_stability_boundary_2d` (`P_c_stat`) → `mirza2025/isotropic_threshold_dimless` (`lambda_crit_0`) at $\beta = 1, \nu = 4$ (exact factor $1/2$ from $A_0 f_A = 2$).
  - *Exact isotropic reductions of Mirza's active-nematic theory*:
    - `rel_mirza2025_disp_eqC7_reduces_ard_eq12_13` (`reduces_to`, err `6.1e-15`): `mirza2025/dispersion_1d_eqC7` (`alpha_plus`) → `kumar2014/ard_linear_stability_eq12_13` (`re_lam_plus`) as $\kappa, \beta \to 0$ with exact 2D $\nu = 4$.
    - `rel_mirza2025_marginal_eqC8_reduces_ard_osc` (`reduces_to`, err `4.9e-09`): `mirza2025/marginal_stability_eqC8` (`lambda_crit`) → `kumar2014/ard_stability_boundary_2d` (`P_c_osc`) as $\kappa \to 0$ ($\nu = 4, \rho = 1$).

---

## Cross-Paper Active-Gel Synthesis Matrix

| Feature | `bois2011` | `kumar2014` | `hannezo2015` | `mietke2019` | `barberi2023` | `mirza2025` |
|---|---|---|---|---|---|---|
| Stress regulators | 1 ($c$) or 2 (ASDM: only $a$ active) | 2 ($A, I$ opposing: $f_A > 0, f_I < 0$) | 1 ($\rho$) | 1 ($c$) + surface geometry | 1 ($C$) + nucleators ($N_a, N_i$) | 1 ($\rho$) + nematic tensor $\mathbf{q}$ |
| Substrate friction $\gamma$ | Finite ($\ell = \sqrt{\eta/\gamma}$) | Finite ($\ell = \sqrt{\eta/\gamma}$) | Finite ($l_h = \sqrt{\eta/\xi}$) | Zero ($\gamma_{\text{ext}} = 0$) | Finite ($\lambda = \sqrt{\eta/\gamma}$) | Finite ($\ell_s = \sqrt{\eta/\gamma}$) |
| 1D/2D longitudinal viscosity | $\eta$ ($\nu = 1$) | $\eta_v$ (1D); $4\eta$ (2D, $\eta_v=3\eta$) | $\eta$ ($\nu = 1$) | $\eta_s + \eta_b$ ($\bar{\nu} = \nu/4$) | $[2(1-1/d)+\tilde{\eta}/\eta]\eta$ | $4\eta$ (2D compressible) |
| Turnover length $\ell_\kappa$ | Only in Eq. 15 / ASDM | $\sqrt{D/\kappa}$ (Eqs. 12–13) | $l_t = \sqrt{D\tau}$ | $\sqrt{D/k}$ | $\sqrt{D_c/k_d}$ | $\ell_D = \sqrt{D/k_d}$ |
| Route to oscillations | None (stationary / ASDM Turing) | Opposing 2-species stress + $\alpha < 1$ or $\rho < 1$ | None (stationary Type-I rings) | Shape–concentration feedback on cylinder | Bistable nucleator chemistry + nonmonotonic $\Pi(C)$ | Isotropic hexagons / nematic bands |
| Empirical `origin: measured` | 0 | 0 | 35 (*Drosophila* trachea) | 0 | 0 | 0 |
