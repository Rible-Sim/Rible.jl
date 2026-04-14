<div align="center">
  <img src="docs/src/assets/logo.svg" height="120" alt="Rible Logo"/>
</div>

<div align="center">
    <h2>Rible.jl</h2>
    <h3>Towards a differentiable simulator for computational design of rigid-flexible robots.</h3>
</div>

<div align="center">

[![Latest Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://docs.rible.dev)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

</div>

---

## Key Features

| Feature | Description |
|:---:|:---|
| ⌨️ **Coordinates** | NCF (Natural Coordinates) - 2D/3D variants with nonminimal formulation. |
| 💥 **Modeling** | Comprehensive library bodies and joints. |
| ⚡ **Dynamics** | Efficient solvers (default: Zhong06 family), contact handling, and linearization, support for (adjoint) sensitivity analysis. |
| 🤖 **Visualization** | Integrated and extensible visualization tools (via Makie). |
| 🛠️ **Extensibility** | Modular design allowing easy extension of new components. |

## Monorepo Structure

Rible.jl is a monorepo structured around a core package and strictly coupled extensions that expand its capabilities.

### 📦 Main Package: Rible.jl

The core package `Rible` provides the fundamental **Multibody Dynamics** engine and essential tools for simulation.

### 🔗 Extensions Architecture

Various extension packages integrate with the core `Rible.jl` ecosystem. 

<div align="center">

| **Package** | **Description** |
| :--- |  :--- |
| 📦 **[Rible.jl](./)**  | The core package for rigid-flexible body dynamics. |
| **Coordinates & Bodies** | | |
| └ 📦 [RibleQCF.jl](./RibleQCF) | Quaternion-based coordinate formulations. |
| **Structures** | | |
| └ 📦 [RibleTensegrity.jl](./RibleTensegrity)  | Tensegrity structures modeling and analysis. |
| **Dynamics Solvers** | | |
| └ 📦 [RibleExtraIntegrators.jl](./RibleExtraIntegrators)  | Include extra integrators that build on the core package and shipped extensions. |


</div>

## Installation

```julia
import Pkg
Pkg.add(url="https://github.com/Rible-Sim/Rible.jl")
```

## Documentation

See the [Getting Started](https://docs.rible.dev/en/get_started/) guide for a walkthrough of your first simulation, or browse the [full documentation](https://docs.rible.dev) for API reference and advanced topics. 

中文用户请参阅[中文文档](https://docs.rible.dev/cn/)。

## Contributing


### Development setup

Clone the repository, enter the checkout, and instantiate the project:

```bash
git clone https://github.com/Rible-Sim/Rible.jl.git
cd Rible.jl
julia --project=. -e 'import Pkg; Pkg.instantiate()'
```

### Validation before opening a pull request

Run the root test suite from the repository root:

```bash
julia --project=. -e 'import Pkg; Pkg.test()'
```

