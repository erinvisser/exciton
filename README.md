# exciton

Two-component exciton model for pre-equilibrium nuclear reactions. Computes internal transition rates (lambda-plus, lambda-conversion) and particle-hole state densities for the Ar40(n,x) reaction, designed to reproduce TALYS output.

## Physics

The exciton model describes the pre-equilibrium phase of a nuclear reaction, where the projectile + target form a compound system that evolves through successive particle-hole (exciton) states before reaching equilibrium.

### Two-component formalism

Protons and neutrons are tracked independently. An exciton state is defined by four integers:

```
(p_pi, h_pi, p_nu, h_nu)  =  (proton particles, proton holes, neutron particles, neutron holes)
```

Total exciton number: `n = p_pi + h_pi + p_nu + h_nu`

### Internal transition rates

The system evolves via:

| Rate | Symbol | Process | Exciton change |
|------|--------|---------|----------------|
| Proton pair creation | `λ⁺_π` | (pπ, hπ) + 1 | n → n + 2 |
| Neutron pair creation | `λ⁺_ν` | (pν, hν) + 1 | n → n + 2 |
| Proton→neutron conversion | `λ_π→ν` | pπ→pν, hπ→hν | n unchanged |
| Neutron→proton conversion | `λ_ν→π` | pν→pπ, hν→hπ | n unchanged |

Each rate is computed via either:
- **Analytical** (closed-form integral, Eq. 13.27 from Koning-Duijvestijn)
- **Numerical** (direct integration over energy, Eqs. 13.17-13.18)

Both analytical and numerical paths include Pauli blocking corrections, finite well depth effects, and pairing corrections (Eqs. 13.5-13.11).

### Matrix element `|M|²`

Two formulations are implemented on different branches:

| Approach | Formula | Branch |
|----------|---------|--------|
| **Matrix element** (simple) | `\|M\|² ∝ A_p/A³ · (7.48 + 4.62e5/(10.7 + E_tot/(n·A_p))³)` | `main`, `MatrixSquared` |
| **Optical model** (OMP) | Above + transmission-coefficient weighting via Koning-Delaroche OMP | `main` |

The `M²` formula includes a tunable 1.20× factor applied only for `PreeqMode::Analytical` (matching TALYS `preeqmode=1`). When n=1 forces analytical mode in numerical mode, this factor is divided back out (`transition_rates.cpp:646-649`).

### State densities

Particle-hole state densities `ω(p_π, h_π, p_ν, h_ν, E_x)` follow the Williams formula with pairwise-equidistant-spacing model, including:
- Pauli blocking corrections
- Finite well depth corrections (Eq. 13.10)
- Superfluid pairing model (Eqs. 13.5-13.7)
- Spin distribution factor (Eqs. 13.59-13.61)

## Branch structure

Three branches document the project's evolution:

```
main  ─── fa96798 (initial repo, full OMP code)
          │
          └── MatrixSquared  ── 854f5ed (strip OMP) ── 01d95bd (re-strip) ── ... ── 96b898f (validated against TALYS)
```

### `main`
Full implemention including the Koning-Delaroche optical model (`CollisionKernel::OpticalModel`), MARLEY MassTable for separation energies, and OMP-weighted collision terms. Intended as a reference implementation matching TALYS' most detailed formalism.

### `MatrixSquared`
OMP stripped, returning to a pure matrix-element-squared kernel (`CollisionKernel::MatrixElement`). Created because:
- The OMP approach introduced fragility — small changes to OMP parameters broke the simple matrix-element agreement
- Maintaining both approaches in the same codebase became unwieldy
- The matrix-element-only version provides a cleaner validation baseline against TALYS

### `OMP`
Originally for OMP-only experiments; currently at the same state as `main`.

**Current work:** The `MatrixSquared` branch has been validated against TALYS output. 81/84 transition rate values agree within <1% relative error. The n=1 discrepancy (originally ~16.7%) was traced to a missing `result /= 1.20` divider dropped during OMP stripping and has been fixed (`transition_rates.cpp:646-649`). Row 0 (n=1) now agrees within 0.014%. The only significant remaining outlier is row 20 (5,5,1,0) lambdanuplus at ~13.4%, understood as a limitation of the simplified `|M|²` formula for extreme exciton partitions.

## Code architecture

```
exciton/
├── include/
│   ├── preeq_common.hh          # Constants, structs, common functions
│   ├── state_density.hh         # Particle-hole state density
│   └── transition_rates.hh      # Rate function declarations
├── src/
│   ├── main_lambda.cpp          # Driver: computes transition rates for Ar40(n,x)
│   ├── main_omega.cpp           # Driver: computes state densities + spin distributions
│   ├── preeq_common.cpp         # g_p, g_n, Pauli correction, finite well, pairing
│   ├── state_density.cpp        # ω(p_π, h_π, p_ν, h_ν, E_x) implementation
│   └── transition_rates.cpp     # λ⁺, λ_conv implementations, M², OMP terms
├── executables/
│   ├── lambda                   # Built binary (transition rates)
│   └── omega                    # Built binary (state densities)
├── docs/
│   └── omp_implementation.html  # Notes on OMP implementation
├── lambda_printouts/            # Validation output for transition rates
│   ├── MatrixSquared/           # Output from MatrixSquared branch
│   ├── OMP/                     # Output from OMP runs (empty)
│   └── *.txt                    # Historical validation runs
├── omega_printouts/             # Validation output for state densities
├── compare_lambda.py            # Compare transition rate output vs reference
├── compare_omega.py             # Compare state density output vs reference
├── Makefile                     # Build system
└── README.md                    # This file
```

### Key equation references (Koning-Duijvestijn 2004)

| Eq. | Description | File |
|-----|-------------|------|
| (13.4) | Particle-hole state density | `state_density.cpp:18` |
| (13.5) | Pairing correction | `preeq_common.cpp:98` |
| (13.6) | Critical exciton number | `preeq_common.cpp:129` |
| (13.7) | Pairing energy gap | `preeq_common.cpp:122` |
| (13.8) | Pauli correction | `preeq_common.cpp:61` |
| (13.9) | Single-particle state density | `preeq_common.cpp:6,12` |
| (13.10) | Finite well depth | `preeq_common.cpp:71` |
| (13.11) | Potential well depth | `preeq_common.cpp:36` |
| (13.17) | Numerical pair creation rate | `transition_rates.cpp:140-352` |
| (13.18) | Numerical conversion rate | `transition_rates.cpp:409-579` |
| (13.24) | M² matrix element (type-weighted) | `transition_rates.cpp:48` |
| (13.26) | M² scaling formula | `transition_rates.cpp:34` |
| (13.27) | Analytical pair creation + conversion | `transition_rates.cpp:86-138, 354-407` |
| (13.28) | Pauli blocking barrier B | `transition_rates.cpp:67` |
| (13.29-30) | Volume-averaged OMP potential | `transition_rates.cpp:12-32` |
| (13.59) | Spin distribution | `preeq_common.cpp:148` |

## Dependencies

### Required
- **C++17** compiler (tested with GCC 13+)
- **MARLEY** (Monte Carlo Archive for Reaction Likelihoods) — private fork
  - `KoningDelarocheOpticalModel` — OMP implementation
  - `Integrator` — Clenshaw-Curtis numerical integration
  - `MassTable` — separation energies for OMP path
  - **Version requirement:** Must be at tag `v1.2.1-exciton` (branch `v121_kd`).
    The `main` branch of MARLEY uses "KDUQ Federal" OMP parameters
    (arXiv:2211.07741) that differ from the original KD03 values used by TALYS 1.6.
    Using the wrong MARLEY version will silently change OMP-weighted transition
    rates and break TALYS agreement on the `main` exciton branch.
- **GSL** (GNU Scientific Library) — linked but currently unused directly
- **ROOT** — for `Math::Integrator` (legacy) and `TGraph` comparison tools

### Build from source

```bash
# Clone MARLEY private fork alongside exciton
git clone <marley-private-url> ../marley-private

# Build
cd exciton
make

# Specific targets
make executables/lambda     # Transition rate calculator
make executables/omega      # State density calculator
make run                    # Run lambda with ARGS
make clean                  # Remove binaries
```

The Makefile expects:
- MARLEY at `../marley-private/` with `.VERSION` file
- ROOT at default prefix (via `root-config`)
- GSL at `/opt/homebrew/lib` (macOS Homebrew) — adjust `LDFLAGS` for your system

### macOS (Homebrew) setup
```bash
brew install gsl root
# MARLEY private fork requires separate setup
```

### Linux setup
Adjust `LDFLAGS` and `ROOT_CFLAGS` in `Makefile` to match your system paths.

## Usage

### Transition rates (`executables/lambda`)

Computes internal transition rates for Ar40(n,x) nuclear reaction with Z=18, N=23 (compound), E_comp=25.6 MeV, V=38 MeV.

```bash
# Numerical mode (default) — Clenshaw-Curtis integration with guard bounds
./executables/lambda

# Analytical mode (closed-form integrals)
./executables/lambda --analytical

# Midpoint integration with 100 bins
./executables/lambda --method midpoint --midbins 100

# No guard bounds
./executables/lambda --no-guards

# Optical model kernel (main branch only)
./executables/lambda --kernel optical-model

# Pass arguments via make
make run ARGS="--analytical --no-guards"
```

**Output format:** 21 exciton states, 8 columns in YANDF-0.4 format:
```
p(p)  h(p)  p(n)  h(n)  λ⁺_π [s⁻¹]  λ⁺_ν [s⁻¹]  λ_π→ν [s⁻¹]  λ_ν→π [s⁻¹]
```

### State densities (`executables/omega`)

Computes particle-hole state densities ω(p_π, h_π, p_ν, h_ν, E_x) for Ar41 and spin distribution factors.

```bash
./executables/omega
```

Output includes state densities for 9 configurations across E_x = 1–25 MeV, and spin distributions for n = 1–12, J = 0–8.

### Comparison scripts

Compare code output against TALYS reference:

```bash
python3 compare_lambda.py \
  lambda_printouts/MatrixSquared/Output/reverted_to_matrixelems.txt \
  /path/to/talys/lambda.out \
  comparison_report.txt \
  1e-4

python3 compare_omega.py \
  omega_printouts/omega_neutron_Ar40_1.txt \
  /path/to/talys/omega.out \
  omega_comparison.txt \
  1e-4
```

## Validation history

The following summarizes key validation runs comparing code output against TALYS reference values. All runs target Ar40(n,x) at E_n = 20 MeV.

### State densities (`executables/omega`)

| Date | Config | Max rel err | Outcome | File |
|------|--------|-------------|---------|------|
| 2026-06-15 | Fixed V=38, no pairing, Clenshaw-Curtis | <1% | PASS | `omega_neutron_Ar40_fixedV_1.txt` |
| 2026-06-15 | With pairing | ~5% at low E_x | Known discrepancy | `omega_neutron_Ar40_1.txt` |
| 2026-06-15 | Spin distribution, talysCompat=true | <0.1% | PASS | `omega_neutron_Ar40_spin*.txt` |
| 2026-06-15 | Without pairing, with guard bounds | <0.5% | PASS | `omega_neutron_Ar40_nopairing_withguard.txt` |

### Transition rates (`executables/lambda`)

| Date | Config | n=1 rel err | Overall | File |
|------|--------|-------------|---------|------|
| 2026-06-16 | Matrix element, analytical, midpoint | ~90% | FAIL | `lambda_comparison_report_analytical.txt` |
| 2026-06-16 | Matrix element, numerical, Clenshaw-Curtis | ~90% | FAIL | `lambda_comparison_report_numerical_clenshaw.txt` |
| 2026-06-22 | Matrix element, numerical, guarded, midBins=50 | 16.6% | Most <1%; 70/84 > 1e-4 abs | `MatrixSquared/Comparison/reverted_to_matrixelems.txt` |
| 2026-06-22 | Fix: `result /= 1.20` for n=1 numerical | **0.014%** | **81/84 < 1%; 3/84 > 1%** | `MatrixSquared/Comparison/reverted_to_matrixelems.txt` |

The n=1 discrepancy (~16.7%) was traced to the `result /= 1.20` divider being accidentally dropped during OMP stripping (commit `854f5ed`). Restored at `transition_rates.cpp:646-649`. After fix, row 0 relative error dropped to 0.014%.

The remaining outliers are: row 20 (5,5,1,0) lambdanuplus ~13.4% — understood as a limitation of the simplified `|M|²` formula for extreme exciton partitions where TALYS' transmission-coefficient approach differs most; and rows 1 and 2 at ~1.7-2.0% — believed to be secondary effects from the same formula difference, not yet investigated in detail.

## Key decisions

### Why separate branches for OMP vs. matrix element?

The OMP approach attempted to match TALYS' bonetti.f90 by weighting collision terms with transmission coefficients from the Koning-Delaroche optical model. While this is more physically complete, it introduced:
- Dependency on MARLEY's OMP implementation (upstream changes could break agreement)
- Fragile calibration (small OMP parameter changes had outsized effects on rates)
- Difficulty isolating whether discrepancies came from the state density, matrix element, or OMP terms

Creating separate branches allows the matrix-element-only approach to be validated cleanly against TALYS, with OMP as a potential future improvement path.

### Why `result /= 1.20` for n=1 in numerical mode?

TALYS' `matrix.f90` applies a 1.20× `|M|²` factor only for `preeqmode=1` (pure analytical mode). The C++ code always routes n=1 through analytical helpers (since numerical integration is degenerate for n=1). These helpers call `M2_element()` without an explicit mode, defaulting to `PreeqMode::Analytical`, which triggers the 1.20× factor. The divider cancels this for numerical-mode runs, matching TALYS behavior.

### Why N=23 (compound) not N=22 (target)?

The initial code used N=22 (the target N for Ar40), but the compound nucleus after neutron absorption is Ar41 with N=23. This was caught during validation; A_compound = Z + N (not Z + N + A_p) was also corrected so `|M|²` uses the compound mass number. See `main_lambda.cpp:38-39` and `transition_rates.cpp:98`.

### Why midBins=50 and guarded Clenshaw-Curtis by default?

The default integration settings (midBins=50 for midpoint, guard bounds for Clenshaw-Curtis) were chosen to give sub-0.5% numerical error, well below the ~1% target agreement with TALYS. Midpoint with midBins=20 is faster but coarser. Guard bounds enforce Pauli-allowed integration limits (Eq. 13.8 minimum energy constraints), preventing unphysical negative contributions.

## References

- Koning, A.J. and Duijvestijn, M.C., "A global pre-equilibrium analysis from 7 to 200 MeV based on the optical model potential", *Nuclear Physics A*, 744:15-76, 2004.
- Kalbach, C., "The Griffin model, complex particles and direct nuclear reactions", *Zeitschrift für Physik A*, 283:401-411, 1977.
- Koning, A.J. and Delaroche, J.P., "Local and global nucleon optical models from 1 keV to 200 MeV", *Nuclear Physics A*, 713:231-310, 2003.
