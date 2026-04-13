using Test
using SafeTestsets

@testset "RibleExtraIntegrators" begin
    @safetestset "Slider Crank" begin include("slider_crank.jl") end
end
