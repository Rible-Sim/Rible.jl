TestEnv.activate()
using Literate
using GLMakie

# Literate.notebook(
#     joinpath(@__DIR__,"pointmass.jl"), 
#     joinpath(@__DIR__,"../../../notebooks/");
#     name="pointmass",
#     execute = true,
#     # documenter = false,
# )

# Literate.notebook(
#     joinpath(@__DIR__,"spinning_top.jl"), 
#     joinpath(@__DIR__,"../../../notebooks/");
#     name="spinning_top",
#     execute = true,
#     # documenter = false,
# )


# Literate.notebook(
#     joinpath(@__DIR__,"woodpecker.jl"), 
#     joinpath(@__DIR__,"../../../notebooks/");
#     name="woodpecker",
#     execute = true,
#     # documenter = false,
# )


## Literate.notebook(
##     joinpath(@__DIR__,"meteor_hammer.jl"), 
##     joinpath(@__DIR__,"../../../notebooks/");
##     name="meteor_hammer",
##     execute = true,
##     # documenter = false,
## )


Literate.notebook(
    joinpath(@__DIR__,"superball.jl"), 
    joinpath(@__DIR__,"../../../notebooks/");
    name="superball",
    execute = true,
    # documenter = false,
)
GLMakie.closeall()