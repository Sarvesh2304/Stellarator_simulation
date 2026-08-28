# Stellarator Physics - Julia Scaffold

A Julia package for self-study calculations involving simplified stellarator
magnetic fields and transport.

## Scope and limitations

This is scaffolding with simplified placeholder physics, not a validated
stellarator code. The transport model does not implement effective helical
ripple (`epsilon_eff`) and has not been benchmarked against VMEC, DKES, SFINCS,
NEO, experiments, or published numerical data. Do not use its numbers for
design or research decisions.

The optimisation objectives are illustrative functions, not validated proxies
for quasi-symmetry, quasi-isodynamicity, or magnetic-well depth. The package is
valuable as a worked exercise in structuring a physics package in Julia.

## Contents

| Module | Intent |
|---|---|
| `MagneticField3D.jl` | 3D field representation and field-line tracing |
| `NeoclassicalTransport.jl` | Simplified transport estimates |
| `StellaratorOptimization.jl` | Illustrative optimisation objectives |
| `TokamakComparison.jl` | Comparison harness |
| `Visualization3D.jl` | Field-line and surface plotting |

## Installation

```bash
git clone https://github.com/Sarvesh2304/Stellarator_simulation.git
cd Stellarator_simulation
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Then run the tests:

```bash
julia --project=. test/test_basic.jl
```

Requires Julia 1.8 or later. Plotting dependencies may require a graphical
backend in headless environments.

## Roadmap

The next research-grade step would be implementing `epsilon_eff` and comparing
the 1/nu neoclassical energy-flux threshold described by Beidler et al.,
*Nature* 596, 221-226 (2021). That work is not included in this scaffold.

## License

MIT. See [LICENSE](LICENSE).

## Author

Sarvesh Pardeshi ([@Sarvesh2304](https://github.com/Sarvesh2304))

Corrections are welcome. Please open an issue with a reproducible example and
the source supporting any proposed physics change.
