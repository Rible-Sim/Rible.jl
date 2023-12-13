TestEnv.activate()
using Literate

Literate.notebook(
    joinpath(@__DIR__,"pointmass.jl"), 
    joinpath(@__DIR__,"../../../notebooks/");
    name="pointmass",
    execute = true,
    # documenter = false,
)
