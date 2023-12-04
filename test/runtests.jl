using Rible
using Test
using SafeTestsets

const GROUP = get(ENV,"GROUP","ALL")

@time begin
    if GROUP == "ALL" || GROUP == "UNIT"
        # unit tests
        @time @safetestset "utils Test" include("unit/loci.jl")
        # @time @safetestset "utils Test" include("unit/pointmass.jl")
        # @time @safetestset "utils Test" include("unit/pointmass.jl")
        # @time @safetestset "utils Test" include("unit/pointmass.jl")
        # @time @safetestset "utils Test" include("unit/pointmass.jl")
        # @time @safetestset "utils Test" include("unit/pointmass.jl")
    end

    if GROUP == "ALL" || GROUP == "MODELING"
        # interface / modeling tests
        @time @safetestset "bodies Test" include("modeling/bodies.jl")
        @time @safetestset "robots Test" include("modeling/robots.jl")
    end

    if GROUP == "ALL" || GROUP == "MECHANICS"
        # mechanics tests
        @time @safetestset "inverse_statics tests"                   include("mechanics/pointmass.jl")
        @time @safetestset "forward_statics tests"                   include("mechanics/pointmass.jl")
        @time @safetestset "stiffness/stability/linearization tests" include("mechanics/pointmass.jl")
        @time @safetestset "adjoint statics dynamics tests"          include("mechanics/pointmass.jl")
        @time @safetestset "dynamic_relax tests"                     include("mechanics/pointmass.jl")
        @time @safetestset "dynamics tests"                          include("mechanics/pointmass.jl")
        @time @safetestset "adjoint dynamics tests"                  include("mechanics/pointmass.jl")
        @time @safetestset "adjoint contact dynamics tests"          include("mechanics/pointmass.jl")
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


