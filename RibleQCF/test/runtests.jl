using Rible
import Rible as RB
using Test
using SafeTestsets

const GROUP = get(ENV,"GROUP","ALL")
if get(ENV, "BENCHMARK", "0") == "1"
    # @info "Running benchmarks (tests skipped)"
    # include("benchmark/run.jl")
else
    @time begin
        if GROUP == "ALL" || GROUP == "UNIT"
            # unit tests
            @time @safetestset "utils Test" include("QCF.jl")
        end

        if GROUP == "ALL" || GROUP == "DYNAMICS"
            # demos tests/ regression tests
            @time @safetestset "Spinning Top Demo Test" include("dynamics.jl")
            # @time @safetestset "Demos Tests" include("demos/nonsmooth/meteor_hammer.jl")
            # @time @safetestset "Demos Tests" include("demos/nonsmooth/pecking_bird.jl")
            # @time @safetestset "Demos Tests" include("demos/nonsmooth/slider_crank.jl")
        end

        if GROUP == "ALL" || GROUP == "ADJOINT"
            # @time @safetestset "adjoint contact dynamics tests"          include("adjoint.jl")
            # @time @safetestset "Policy Gradients Test" include("policy_gradients_test.jl")
        end

    end
end
