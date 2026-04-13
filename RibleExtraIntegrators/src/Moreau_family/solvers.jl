
struct Moreau{T} <: AbstractIntegrator
    θ::T
end

include("constant_mass.jl")
include("ccp_constant_mass.jl")
