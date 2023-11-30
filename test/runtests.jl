using Rible
using Test
using SafeTestsets

@testset "Rible.jl" begin
    # Write your tests here.
end


@safetestset "utils Test" begin include("demos/nonsmooth/pointmass.jl") end


# demos
@safetestset "Demos Tests" begin include("demos/nonsmooth/pointmass.jl") end
@safetestset "Demos Tests" begin include("demos/nonsmooth/spinning_top.jl") end
@safetestset "Demos Tests" begin include("demos/nonsmooth/meteor_hammer.jl") end
@safetestset "Demos Tests" begin include("demos/nonsmooth/superball.jl") end
@safetestset "Demos Tests" begin include("demos/nonsmooth/pecking_bird.jl") end
@safetestset "Demos Tests" begin include("demos/nonsmooth/slider_crank.jl") end
