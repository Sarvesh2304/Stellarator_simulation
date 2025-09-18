
# Stellarator Physics Analysis Package

A comprehensive Julia package for stellarator fusion reactor physics analysis, featuring 3D magnetic field calculations, neoclassical transport modeling, optimization algorithms, and comparison with tokamak performance.

## Features

### 🔬 Core Physics
- **3D Magnetic Field Calculations**: Complete stellarator magnetic field modeling with Fourier harmonics
- **Neoclassical Transport**: Stellarator-specific transport coefficient calculations
- **Boozer Coordinates**: Advanced coordinate system for stellarator analysis
- **Magnetic Surface Analysis**: Calculation of magnetic surfaces and safety factors

### ⚡ Optimisation Algorithms
- **Quasi-Symmetry Optimisation**: Minimize magnetic field asymmetry
- **Quasi-Isodynamicity Optimisation**: Optimise for isodynamic magnetic fields
- **Magnetic Well Optimisation**: Maximize magnetic well depth for stability
- **Transport Optimisation**: Minimize neoclassical transport

### 🔄 Tokamak Comparison
- **Performance Metrics**: Direct comparison of stellarator vs tokamak performance
- **Transport Analysis**: Side-by-side transport coefficient comparison
- **Stability Analysis**: Comparison of stability properties and beta limits
- **Confinement Scaling**: Analysis of confinement time scaling laws

### 📊 3D Visualisation
- **Magnetic Field Lines**: Interactive 3D field line tracing
- **Plasma Surfaces**: 3D visualisation of magnetic surfaces
- **Transport Profiles**: 2D and 3D transport coefficient visualisation
- **Optimisation Results**: Real-time optimisation progress visualisation

## Installation

1. Clone the repository:
```bash
git clone <repository-url>
cd "Stellarator Physics"
```

2. Start Julia and activate the environment:
```julia
julia> using Pkg
julia> Pkg.activate(".")
julia> Pkg.instantiate()
```

3. Load the package:
```julia
julia> using StellaratorPhysics
```

## Quick Start

```julia
using StellaratorPhysics

# Create a stellarator configuration
R₀ = 1.0  # Major radius [m]
a = 0.2   # Minor radius [m]
N = 5     # Number of field periods
B₀ = 1.0  # Reference magnetic field [T]

stellarator_bfield = MagneticField3D(R₀, a, N, B₀)

# Set up plasma parameters
T_e = 1000.0  # Electron temperature [eV]
T_i = 1000.0  # Ion temperature [eV]
n_e = 1e20    # Electron density [m^-3]
n_i = 1e20    # Ion density [m^-3]

transport = NeoclassicalTransport(stellarator_bfield, T_e, T_i, n_e, n_i)

# Calculate transport coefficients
s = 0.5  # Normalized radius
coeffs = calculate_transport_coefficients(transport, s)

# Create 3D visualization
plot_3d = plot_plasma_surfaces(stellarator_bfield, [0.2, 0.4, 0.6, 0.8])
```

## Examples

### Complete Analysis Example
Run the comprehensive example:
```julia
include("examples/stellarator_analysis_example.jl")
```

This example demonstrates:
- 3D magnetic field calculations
- Neoclassical transport analysis
- Stellarator optimization
- Comparison with tokamak performance
- 3D visualization

### Individual Module Examples

#### Magnetic Field Analysis
```julia
# Calculate magnetic field at a point
R, φ, Z = 1.2, 0.0, 0.1
B_R, B_φ, B_Z = calculate_magnetic_field(stellarator_bfield, R, φ, Z)

# Trace a field line
field_line = trace_field_line(stellarator_bfield, R, φ, Z, 10.0)

# Find magnetic surface
surface = find_magnetic_surface(stellarator_bfield, 0.5)
```

#### Transport Analysis
```julia
# Calculate transport coefficients
coeffs = calculate_transport_coefficients(transport, 0.5)

# Calculate bootstrap current
j_bs = bootstrap_current(transport, 0.5)

# Calculate radial electric field
E_r = radial_electric_field(transport, 0.5)
```

#### Optimisation
```julia
# Optimise for quasi-symmetry
result, optimal_harmonics = optimize_quasi_symmetry(stellarator_bfield, 1000)

# Multi-objective optimization
weights = Dict(
    "quasi_symmetry" => 0.4,
    "magnetic_well" => 0.3,
    "transport" => 0.3
)
result = multi_objective_optimization(stellarator_bfield, weights, 1000)
```

#### Comparison with Tokamak
```julia
# Create tokamak field
tokamak_bfield = create_tokamak_field(R₀, a, 1.0, 2.0, B₀)
tokamak_transport = NeoclassicalTransport(tokamak_bfield, T_e, T_i, n_e, n_i)

# Compare performance
comparison = TokamakComparison(stellarator_bfield, tokamak_bfield, 
                              stellarator_transport, tokamak_transport)
results = comprehensive_comparison(comparison, 0.5)
```

## Package Structure

```
StellaratorPhysics/
├── src/
│   ├── StellaratorPhysics.jl      # Main module
│   ├── MagneticField3D.jl         # 3D magnetic field calculations
│   ├── NeoclassicalTransport.jl   # Transport modeling
│   ├── StellaratorOptimization.jl # Optimisation algorithms
│   ├── TokamakComparison.jl       # Tokamak comparison
│   ├── Visualization3D.jl         # 3D visualization
│   └── utils.jl                   # Utility functions
├── examples/
│   └── stellarator_analysis_example.jl
├── Project.toml
└── README.md
```

## Dependencies

- **Julia 1.8+**: Required for optimal performance
- **Plots.jl**: For 2D plotting and visualization
- **PlotlyJS.jl**: For interactive 3D visualization
- **Optim.jl**: For optimization algorithms
- **NLopt.jl**: For advanced optimization
- **DifferentialEquations.jl**: For field line tracing
- **ForwardDiff.jl**: For automatic differentiation
- **SpecialFunctions.jl**: For special mathematical functions

## Physics Background

### Stellarator vs Tokamak
Stellarators are toroidal fusion devices that use 3D magnetic field configurations to confine plasma, unlike tokamaks which rely on axisymmetric fields and plasma current. Key advantages include:
- **Steady-state operation**: No need for plasma current drive
- **Reduced MHD instabilities**: 3D field provides additional stability
- **Flexible design**: Can optimize for various physics objectives

### Neoclassical Transport
Neoclassical transport in stellarators differs from tokamaks due to:
- **3D magnetic field geometry**: Affects particle orbits and transport
- **Magnetic field ripple**: Creates additional transport channels
- **Quasi-symmetry**: Can reduce transport to tokamak-like levels

### Optimisation Objectives
- **Quasi-Symmetry**: Minimize magnetic field asymmetry to reduce transport
- **Quasi-Isodynamicity**: Optimise for isodynamic magnetic fields
- **Magnetic Well**: Maximize magnetic well depth for stability
- **Transport**: Minimize neoclassical transport coefficients

## Contributing

Contributions are welcome! Please see the contributing guidelines for details on:
- Code style and formatting
- Testing requirements
- Documentation standards
- Pull request process

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this package in your research, please cite:

```bibtex
@software{stellarator_physics,
  title={Stellarator Physics Analysis Package},
  author={Stellarator Research Team},
  year={2024},
  url={https://github.com/your-repo/stellarator-physics}
}
```

## Acknowledgments

- Based on stellarator physics theory and VMEC code
- Inspired by the work of the stellarator research community
- Built with the Julia scientific computing ecosystem

## Support

For questions, issues, or contributions:
- Open an issue on GitHub
- Contact the development team
- Check the documentation and examples

---

**Note**: This package is designed for research and educational purposes. For production fusion reactor design, consult with fusion physics experts and use validated codes.
