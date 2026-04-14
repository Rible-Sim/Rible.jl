module RibleRungeKuttaExt

import Rible as RB
import RungeKutta as RK
using ProgressMeter
using Printf
using LinearAlgebra

import Rible: AbstractIntegrator, AbstractDynamicsProblem, Simulator, DynamicsSolver, AbstractCoordinatesState, CoordinatesState
import Rible: build_mass_matrices, interpolate!, guess_newton!, populate!, interpolate!, get_num_of_actions, get_numbertype
import Rible: gen_force!,  gen_force_state_jacobian!, cstr_forces_jacobian!, cstr_jacobian!, cstr_function!
import Rible: generate_cache, solve!
import Rible: compute_rungekutta_residual!, compute_rungekutta_jacobian!

import Rible: RKIntegrator, RungeKuttaConstants, RungeKuttaJacobianWorkspace, RungeKutta_Constant_Mass_Cache, RungeKuttaSolverState



include("primal_solver.jl")

end