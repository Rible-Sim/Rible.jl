
struct GeneralizedAlpha{T} <: AbstractIntegrator
    αm::T
    αf::T
    γ::T
    β::T
end


function GeneralizedAlpha(ρ∞)
    αm = (2ρ∞-1)/(ρ∞+1)
    αf = ρ∞/(ρ∞+1)
    γ = 1/2 + αf - αm
    β = 1/4*(γ+1/2)^2
    GeneralizedAlpha(αm,αf,γ,β)
end

function GeneralizedAlpha(ρ∞,h)
    αm = (2ρ∞-1)/(ρ∞+1)
    αf = ρ∞/(ρ∞+1)
    γ = 1/2 + αf - αm
    β = 1/4*(γ+1/2)^2
    γₜ = (1-αm)/(1-αf)/(γ*h)
    βₜ = h*β/γ - h/2
    GeneralizedAlpha(αm,αf,γ,β,γₜ,βₜ)
end

function Newmark()
    αf = αm = 0.0
    γ = 1/2
    β = 1/4
    GeneralizedAlpha(αm,αf,γ,β)
end

include("primal_solver.jl")
include("ccp_solver.jl")

# include("nonsmooth.jl")
# include("nssfc.jl")
