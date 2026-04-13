using Rible
import Rible as RB
using Test
using SafeTestsets

const GROUP = get(ENV,"GROUP","ALL")
const ACTIVE_ROOT_GROUPS = ("UNIT", "DYNAMICS_A", "DYNAMICS_B", "ADJOINT_A", "ADJOINT_B")

group_enabled(groups::AbstractString...) = GROUP == "ALL" || GROUP in groups

if get(ENV, "BENCHMARK", "0") == "1"
    @info "Running benchmarks (tests skipped)"
    include("benchmark/run.jl")
else
    @time begin
        if group_enabled("UNIT")
            # unit tests
            @time @safetestset "utils Test" include("unit/loci.jl")
            @time @safetestset "assets Test" include("unit/assets.jl")
        end

        if group_enabled("MODELING")
            @info "MODELING has no standalone root safetestsets in the current CI baseline; root verification is partitioned into $(join(ACTIVE_ROOT_GROUPS, ", ")), while modeling coverage comes from the active dynamics and package-level suites."
        end

        if group_enabled("DYNAMICS", "DYNAMICS_A")
            # demos tests/ regression tests
            @time @safetestset "Pointmass Demo Test" include("dynamics/pointmass.jl")
        end

        if group_enabled("DYNAMICS", "DYNAMICS_B")
            @time @safetestset "Spinning Top Demo Test" include("dynamics/spinning_top.jl")
            # @time @safetestset "Demos Tests" include("demos/nonsmooth/meteor_hammer.jl")
            # @time @safetestset "Demos Tests" include("demos/nonsmooth/pecking_bird.jl")
            # @time @safetestset "Demos Tests" include("demos/nonsmooth/slider_crank.jl")
        end

        if group_enabled("ADJOINT", "ADJOINT_A")
            @info "Legacy mechanics pointmass safetestsets are intentionally retired from the root CI baseline; active adjoint verification lives in the focused adjoint suites below."
            @time @safetestset "BSpline Basis" include("adjoint/bspline_test.jl")
            @time @safetestset "ApproxFun Policies" include("adjoint/approxfun_policy_test.jl")
        end

        if group_enabled("ADJOINT", "ADJOINT_B")
            @time @safetestset "Policy Gradients" include("adjoint/policy_gradients_test.jl")
            @time @safetestset "adjoint contact dynamics tests"          include("adjoint.jl")
        end

    end
end
