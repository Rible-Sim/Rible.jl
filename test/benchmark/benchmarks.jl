using BenchmarkTools
using Rible
import Rible as RB
using LinearAlgebra
using Logging
using Rotations
using StaticArrays

include(joinpath(@__DIR__, "..", "..", "examples", "robots", "pointmass.jl"))

# Keep the suite lightweight for CI while still catching allocation regressions.
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 0.01
BenchmarkTools.DEFAULT_PARAMETERS.samples = 1
BenchmarkTools.DEFAULT_PARAMETERS.evals = 1
global_logger(SimpleLogger(devnull, Logging.Error))

const POINTMASS_PARAMS = let
    surface = RB.StaticContactSurfaces([RB.HalfSpace([0, 0, 1.0], zeros(3))])
    policy = RB.NoPolicy()
    initial_state = SVector(0.0, 0.0, 1.0, 0.5, 0.0, 0.0)
    (; surface, policy, initial_state)
end

function fresh_pointmass_robot()
    new_pointmass(
        ; origin_position = POINTMASS_PARAMS.initial_state[1:3],
          origin_velocity = POINTMASS_PARAMS.initial_state[4:6],
          term_position = SVector(0.5, 0.0, 0.0),
          μ = 0.1,
          e = 0.5,
    )
end

function make_cost_fixture()
    robot = fresh_pointmass_robot()
    (; robot)
end

const SUITE = BenchmarkGroup()

SUITE["pointmass"] = BenchmarkGroup()
SUITE["pointmass"]["assemble_M"] = @benchmarkable RB.assemble_M(structure) setup=(fixture = make_cost_fixture(); structure = fixture.robot.structure)
