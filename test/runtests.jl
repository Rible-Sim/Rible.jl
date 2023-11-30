using Rible
using Test
using SafeTestsets

const GROUP = get(ENV,"GROUP","ALL")

@time begin
    if GROUP == "ALL" || GROUP == "UNIT"
        # unit tests
        @time @safetestset "utils Test" include("unit/loci.jl")
        # @time @safetestset "utils Test" include("unit/nonsmooth/pointmass.jl")
        # @time @safetestset "utils Test" include("unit/nonsmooth/pointmass.jl")
        # @time @safetestset "utils Test" include("unit/nonsmooth/pointmass.jl")
        # @time @safetestset "utils Test" include("unit/nonsmooth/pointmass.jl")
        # @time @safetestset "utils Test" include("unit/nonsmooth/pointmass.jl")
    end

    if GROUP == "ALL" || GROUP == "MODELING"
        # interface / modeling tests
        @time @safetestset "utils Test" include("modeling/nonsmooth/pointmass.jl")
        @time @safetestset "utils Test" include("modeling/nonsmooth/pointmass.jl")
        @time @safetestset "utils Test" include("modeling/nonsmooth/pointmass.jl")
        @time @safetestset "utils Test" include("modeling/nonsmooth/pointmass.jl")
        @time @safetestset "utils Test" include("modeling/nonsmooth/pointmass.jl")
        @time @safetestset "utils Test" include("modeling/nonsmooth/pointmass.jl")
    end

    if GROUP == "ALL" || GROUP == "MECHANICS"
        # mechanics tests
        # inverse_statics tests
        # forward_statics tests
        # stiffness/stability/linearization tests
        # adjoint statics dynamics tests
        # dynamic_relax tests
        # dynamics tests
        # adjoint dynamics tests
        # contact dynamics tests
        # adjoint contact dynamics tests
        @time @safetestset "utils Test" include("mechanics/nonsmooth/pointmass.jl")
        @time @safetestset "utils Test" include("mechanics/nonsmooth/pointmass.jl")
        @time @safetestset "utils Test" include("mechanics/nonsmooth/pointmass.jl")
        @time @safetestset "utils Test" include("mechanics/nonsmooth/pointmass.jl")
        @time @safetestset "utils Test" include("mechanics/nonsmooth/pointmass.jl")
        @time @safetestset "utils Test" include("mechanics/nonsmooth/pointmass.jl")
    end

    if GROUP == "ALL" || GROUP == "DEMOS"
        # demos tests/ regression tests
        @time @safetestset "Demos Tests" include("demos/nonsmooth/pointmass.jl")
        @time @safetestset "Demos Tests" include("demos/nonsmooth/spinning_top.jl")
        @time @safetestset "Demos Tests" include("demos/nonsmooth/meteor_hammer.jl")
        @time @safetestset "Demos Tests" include("demos/nonsmooth/superball.jl")
        @time @safetestset "Demos Tests" include("demos/nonsmooth/pecking_bird.jl")
        @time @safetestset "Demos Tests" include("demos/nonsmooth/slider_crank.jl")
    end
end


