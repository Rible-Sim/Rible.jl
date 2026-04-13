```@meta
CurrentModule = Rible
ShareDefaultModule = true
```

# API Reference

This page lists the main APIs of `Rible` by functional module, auto-generated from docstrings.

## Coordinates

```@autodocs
Modules = [Rible, Rible.NCF]
Pages = [
    "coordinates/coordinate_interfaces.jl",
    "coordinates/nonminimal_coordinates.jl",
    "coordinates/natural_coordinates/types.jl",
    "coordinates/natural_coordinates/interface_implementations.jl",
    "coordinates/natural_coordinates/constructors.jl",
    "coordinates/natural_coordinates/functions.jl",
    "coordinates/natural_coordinates/constraints.jl",
    "coordinates/natural_coordinates/joints.jl",
    "coordinates/natural_coordinates/mass_matrix.jl",
    "coordinates/natural_coordinates/LNC.jl",
    "coordinates/natural_coordinates/NCF.jl",
    "coordinates/presfree_coordinates.jl",
    "coordinates/contact_state.jl",
    "utils/signifier.jl",
    "utils/utils.jl",
    "utils/rotations.jl",
]
Private = true
```

## Bodies

```@autodocs
Modules = [Rible]
Pages = [
    "abstract_types.jl",
    "body/types.jl",
    "body/body_interfaces.jl",
    "body/coordinate_specific.jl",
    "body/implementations.jl",
]
Private = true
```

## Structure, Constraints, and Connectivity

```@autodocs
Modules = [Rible]
Pages = [
    "structures/interfaces.jl",
    "structures/structure.jl",
    "structures/connectivity.jl",
    "structures/constraints.jl",
    "structures/mass_matrices.jl",
    "structures/linearization.jl",
    "structures/auxiliary.jl",
    "structures/mutate.jl",
    "structures/pres_free_state.jl",
]
Private = true
```

## Force Elements, Joints, and Generalized Forces

```@autodocs
Modules = [Rible]
Pages = [
    "apparatus/apparatus.jl",
    "force/spring_dampers.jl",
    "force/rheonomic_joint.jl",
    "joint/interfaces.jl",
    "joint/joints.jl",
    "mechanics/generalized_forces.jl",
    "mechanics/materials.jl",
    "mechanics/linearized_mechanics.jl",
]
Private = true
```

## Environment and Contact

```@autodocs
Modules = [Rible]
Pages = [
    "environments/env.jl",
    "environments/geometries/rigids.jl",
    "environments/geometries/surfaces.jl",
    "mechanics/contact.jl",
]
Private = true
```

## Control — Policies, Actuators, Gauges, Costs

```@autodocs
Modules = [Rible]
Pages = [
    "policy/basics.jl",
    "policy/basis_openloop.jl",
    "policy/feedback_static_linear.jl",
    "policy/feedback_proportional_derivative.jl",
    "policy/discrete_linear.jl",
    "policy/discrete_pd.jl",
    "control/control.jl",
    "control/reg_actuator.jl",
    "control/policy.jl",
    "control/actuators.jl",
    "control/gauges.jl",
    "control/capta.jl",
    "control/cost.jl",
    "control/cost_helpers.jl",
]
Private = true
```

## Dynamics Solvers and Sensitivity Analysis

```@autodocs
Modules = [Rible]
Pages = [
    "dynamics_solvers/solvers.jl",
    "dynamics_solvers/workspaces.jl",
    "dynamics_solvers/complementarity_solvers.jl",
    "dynamics_solvers/zhong06_res_and_jac.jl",
    "dynamics_solvers/RungeKutta_family/solvers.jl",
    "dynamics_solvers/RungeKutta_family/primal_solver.jl",
    "dynamics_solvers/Zhong06_family/solvers.jl",
    "dynamics_solvers/Zhong06_family/res_and_jac.jl",
    "dynamics_solvers/Zhong06_family/base/primal_solver.jl",
    "dynamics_solvers/Zhong06_family/base/adjoint_solver.jl",
    "dynamics_solvers/Zhong06_family/base/direct_solver.jl",
    "dynamics_solvers/Zhong06_Constant_Mass_family/CCP_mono/primal_solver.jl",
    "dynamics_solvers/Zhong06_Constant_Mass_family/CCP_mono/adjoint_solver.jl",
    "dynamics_solvers/Zhong06_Constant_Mass_family/CCP_mono/direct_solver.jl",
    "dynamics_solvers/Zhong06_Constant_Mass_family/CCP_inner/primal_solver.jl",
]
Private = true
```

## Robot and Post-processing

```@autodocs
Modules = [Rible]
Pages = [
    "robots/robot.jl",
    "postprocess/analysis.jl",
    "postprocess/vis/plot.jl",
    "postprocess/vis/plot_traj.jl",
    "postprocess/vis/mesh.jl",
    "postprocess/vis/recipe.jl",
    "postprocess/vis/theme.jl",
]
Private = true
```

## Visualization Recipes

```@docs
Rible.Vis
Rible.vis
Rible.vis!
```
