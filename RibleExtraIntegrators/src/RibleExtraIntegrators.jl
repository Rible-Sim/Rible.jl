
module RibleExtraIntegrators

using Rible
using LinearAlgebra
using Printf
using NLsolve
using ProgressMeter
using EponymTuples
using StructArrays
using StaticArrays
using FiniteDiff
using FiniteDifferences

import Rible:
    AbstractIntegrator,
    AbstractDynamicsProblem,
    AbstractBodySolver,
    AbstractApparatusSolver,
    AbstractContactSolver,
    DynamicsProblem,
    DynamicsSolver,
    DiscreteAdjointDynamicsSolver,
    Simulator,
    InnerLayerContactSolver,
    MonolithicContactSolver,
    RestitutionFrictionCombined,
    NewtonRestitution,
    CoulombFriction,
    Frictionless,
    StaticContactSurfaces,
    Plane,
    build_mass_matrices,
    assemble_M, assemble_M!, assemble_M⁻¹!, assemble_M⁻¹, assemble_∂Mq̇∂q!, assemble_∂Mq̇∂q, 
    make_cstr_jacobian, cstr_jacobian!,
    make_cstr_function, cstr_function!,
    cstr_forces_jacobian, cstr_forces_jacobian!, 
    cstr_velocity_jacobian!,
    gen_force!, update_bodies!,
    gen_force_state_jacobian!,
    has_constant_mass_matrix,
    get_numbertype, get_num_of_actions, get_num_of_aux_var

import Rible: solve!, interpolate!, populate!, guess_newton!, control!, generate_cache, create_line_search_functions, NewtonWorkspace, Zhong06, Zhong06Constants, Zhong06SolverState, Zhong06JacobianWorkspace

include("Alpha_family/solvers.jl")
include("Moreau_family/solvers.jl")


end
