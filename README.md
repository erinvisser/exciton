# exciton

Two-component exciton model for pre-equilibrium nuclear reactions. Computes internal transition rates (lambda-plus, lambda-conversion) and particle-hole state densities for the Ar40(n,x) reaction, designed to reproduce TALYS output.

**Documentation:** [https://erinvisser.github.io/exciton/](https://erinvisser.github.io/exciton/)

## Overview

The exciton model describes the pre-equilibrium phase of a nuclear reaction, where the projectile + target form a compound system that evolves through successive particle-hole (exciton) states before reaching equilibrium.

This code implements:
- **Transition rates** — proton/neutron pair creation and conversion rates (λ⁺_π, λ⁺_ν, λ_π→ν, λ_ν→π)
- **State densities** — particle-hole state densities ω(p_π, h_π, p_ν, h_ν, E_x) with Pauli corrections, finite well depth, and pairing
- **Two matrix element formulations** — simple |M|² formula (validated) and Koning-Delaroche OMP-weighted kernel

## Branches

| Branch | Status | Description |
|--------|--------|-------------|
| `main` | Full codebase | All branches merged, OMP + matrix element |
| `MatrixSquared` | **Validated** (81/84 < 1%) | Pure matrix-element kernel, stripped of OMP |
| `OMPfromScratch` | In progress | OMP kernel rebuilt from validated M² codebase |
| `OMP` | Reference | Original OMP implementation experiments |

## Quick start

```bash
# Clone with MARLEY dependency
git clone --recurse-submodules <url>
cd exciton
make

# Run transition rates
./executables/lambda

# Run state densities
./executables/omega
```

Requires C++17, MARLEY (v1.2.1-exciton tag), GSL, and ROOT. See the [documentation](https://erinvisser.github.io/exciton/) for full build instructions, physics details, and validation history.
