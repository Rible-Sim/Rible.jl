module ANCF
using LinearAlgebra
using StaticArrays
using Parameters
using ForwardDiff
using DocStringExtensions
import FastGaussQuadrature as FGQ
using QuadGK

"""
所有Absolute Nodal Coordinates abstract type。
"""
abstract type ANC{T} end
"""
所有 two-dimensionalAbsolute Nodal Coordinates abstract type。
"""
abstract type ANC2D{T} <: ANC{T} end
"""
所有 three-dimensionalAbsolute Nodal Coordinates abstract type。
"""
abstract type ANC3D{T} <: ANC{T} end
"""
所有坐标数目为4的 two-dimensionalAbsolute Nodal Coordinates abstract type，用于刚性直杆。
"""
abstract type ANC2D4C{T} <: ANC2D{T} end
"""
所有坐标数目为6的 two-dimensionalAbsolute Nodal Coordinates abstract type，用于任意形状刚体。
"""
abstract type ANC2D6C{T} <: ANC2D{T} end
"""
所有坐标数目为6的 three-dimensionalAbsolute Nodal Coordinates abstract type，用于刚性直杆。
"""
abstract type ANC3D6C{T} <: ANC3D{T} end
"""
所有坐标数目为12的 three-dimensionalAbsolute Nodal Coordinates abstract type，用于 three-dimensional绳索（Low Order Cable）。
"""
abstract type ANC3D12C{T} <: ANC3D{T} end


make_I(T,N) = SMatrix{N,N}(one(T)*I)
const I2 = make_I(Bool,2)
const I3 = make_I(Bool,3)

"""
返回空间维数。
$(TYPEDSIGNATURES)
"""
get_ndim(::ANC2D) = 2
get_ndim(::ANC3D) = 3
"""
返回Absolute Nodal Coordinates所构成的坐标系的维数。
$(TYPEDSIGNATURES)
"""
get_nlocaldim(::ANC2D4C) = 1
get_nlocaldim(::ANC2D6C) = 2
get_nlocaldim(::ANC3D6C) = 1
get_nlocaldim(::ANC3D12C) = 1
"""
返回坐标个数。
$(TYPEDSIGNATURES)
"""
get_ncoords(::ANC2D4C) = 4
get_ncoords(::ANC2D6C) = 6
get_ncoords(::ANC3D6C) = 6
get_ncoords(::ANC3D12C) = 12
"""
返回约束方程个数。
$(TYPEDSIGNATURES)
"""
get_nconstraints(::ANC2D4C) = 0
get_nconstraints(::ANC2D6C) = 0
get_nconstraints(::ANC3D6C) = 0
get_nconstraints(::ANC3D12C) = 0
"""
返回自由度数。
$(TYPEDSIGNATURES)
"""
get_ndof(ancs::ANC) =  get_ncoords(ancs) - get_nconstraints(ancs)


"""
坐标数目为12的 three-dimensionalAbsolute Nodal Coordinates类，使用1个基本点、3个基本向量。
$(TYPEDEF)
"""
struct ANC3DRURU{T,ρT} <: ANC3D12C{T}
    radius::T
    E::T
    I::T
    A::T
    L::T
    ρ::ρT
end
# radius

function ANC3DRURU(ρ::Number;E,radius,L)
    ANC3DRURU((x)->ρ;E,radius,L)
end

function ANC3DRURU(ρ::Function;E,radius,L)
    A = π*radius^2
    I = π*radius^4/4
    ANC3DRURU(radius,E,I,A,L,ρ)
end

function make_x(::ANC3DRURU)
    function x(r̄)
        r̄[begin]
    end
end


# Shape functions
"""
Return shape functions
$(TYPEDSIGNATURES)
"""
function find_ξ(ancs::ANC3DRURU)
    (;L) = ancs
    function ξ(x)
        2*x/L - 1
    end
end

function find_x(ancs::ANC3DRURU)
    (;L) = ancs
    function x(ξ)
        L/2*(ξ + 1)
    end
end

function make_S(ancs::ANC3DRURU)
    (;L) = ancs
    _ξ = find_ξ(ancs)
    function S(x)
        ξ = _ξ(x)
        S_raw = SA[1//2 - 3//4*ξ + 1//4*ξ^3  1//8*L*( 1 - ξ - ξ^2 + ξ^3) 1//2 + 3//4*ξ - 1//4*ξ^3  1//8*L*(-1 - ξ + ξ^2 + ξ^3);]
        kron(S_raw,I3)
    end
end

function make_Sₓ(ancs::ANC3DRURU)
    (;L) = ancs
    _ξ = find_ξ(ancs)
    function Sₓ(x)
        ξ = _ξ(x)
        Sₓ_raw = 2/L.*SA[-3//4 + 3//4*ξ^2  1//8*L*(-1 - 2ξ + 3ξ^2) 3//4 - 3//4*ξ^2  1//8*L*(-1 + 2ξ + 3ξ^2);]
        kron(Sₓ_raw,I3)
    end
end

function make_Sₓₓ(ancs::ANC3DRURU)
    (;L) = ancs
    _ξ = find_ξ(ancs)
    function Sₓₓ(x)
        ξ = _ξ(x)
        Sₓₓ_raw = 4/L^2 .*SA[3//2*ξ  1//8*L*(-2 + 6ξ) -3//2*ξ  1//8*L*(2 + 6ξ);]
        kron(Sₓₓ_raw,I3)
    end
end

# function find_ξ(ancs::ANC3DRURU)
#     (;L) = ancs
#     function ξ(x̄)
#         x̄/L
#     end
# end
function make_r(ancs::ANC3DRURU,e)
    _S = make_S(ancs)
    function r(x)
        _S(x)*e
    end
end

function make_rₓ(ancs::ANC3DRURU,e)
    _Sₓ = make_Sₓ(ancs)
    function rₓ(x)
        _Sₓ(x)*e
    end
end

function make_rₓₓ(ancs::ANC3DRURU,e)
    _Sₓₓ = make_Sₓₓ(ancs)
    function rₓₓ(x)
        _Sₓₓ(x)*e
    end
end

function make_εₓₓ(ancs::ANC3DRURU,e)
    _rₓ = make_rₓ(ancs,e)
    function εₓₓ(x)
        rₓ = _rₓ(x)
        (transpose(rₓ)*rₓ - 1)/2
    end
end

function make_κ(ancs::ANC3DRURU,e)
    _rₓ  = make_rₓ(ancs,e)
    _rₓₓ = make_rₓₓ(ancs,e)
    function κ(x)
        rₓ = _rₓ(x)
        rₓₓ = _rₓₓ(x)
        norm(rₓ×rₓₓ)/norm(rₓ)^3
    end
end

function LinearAlgebra.cross(
        a::StaticArray{Tuple{3}, T1, 1},
        b::StaticArray{Tuple{3,N}, T2, 2}) where {T1,T2,N}
    reduce(hcat,Ref(a).×SVector{N}(eachcol(b)))
end

function make_∂κ∂eᵀ(ancs::ANC3DRURU,e)
    _rₓ  = make_rₓ(ancs,e)
    _rₓₓ = make_rₓₓ(ancs,e)
    _Sₓ  = make_Sₓ(ancs)
    _Sₓₓ = make_Sₓₓ(ancs)
    function ∂κ∂eᵀ(x)
        rₓ  = _rₓ(x)
        rₓₓ = _rₓₓ(x)
        Sₓ  = _Sₓ(x)
        Sₓₓ = _Sₓₓ(x)
        rrₓₓₓ = rₓ × rₓₓ
        f = norm(rrₓₓₓ)
        normrₓ = norm(rₓ)
        g = normrₓ^3
        ∂g∂eᵀ = 3*normrₓ*transpose(Sₓ)*rₓ
        ∂f∂eᵀ = 1/f*transpose(-rₓₓ×Sₓ .+ rₓ×Sₓₓ)*(rrₓₓₓ)
        @. 1/g^2*(g*∂f∂eᵀ-f*∂g∂eᵀ)
    end
end


function make_Q(ancs::ANC3DRURU)
    (;E,A,I) = ancs
    _x = find_x(ancs)
    g3ξ,g3w = FGQ.gausslegendre(3)
    g5ξ,g5w = FGQ.gausslegendre(5)
    g3x = _x.(g3ξ)
    g5x = _x.(g5ξ)
    _Sₓ  = make_Sₓ(ancs)
    _Sₓₓ = make_Sₓₓ(ancs)
    function Q(e)
        _rₓ  = make_rₓ(ancs,e)
        _rₓₓ = make_rₓₓ(ancs,e)   
        # _κ = make_κ(ancs,e)
        # _∂κ∂eᵀ = make_∂κ∂eᵀ(ancs,e)
        # _εₓₓ = make_εₓₓ(ancs,e)
        function qA(x)
            rₓ = _rₓ(x)
            E*A*(transpose(rₓ)*rₓ - 1)/2*transpose(_Sₓ(x))*rₓ
        end
        function qI(x) 
            rₓ = _rₓ(x)
            rₓₓ = _rₓₓ(x)
            normrₓ = norm(rₓ)
            rrₓₓₓ = rₓ × rₓₓ
            f = norm(rrₓₓₓ)
            g = normrₓ^3
            κ = f/g
            Sₓ  = _Sₓ(x)
            Sₓₓ = _Sₓₓ(x)
            ∂g∂eᵀ = 3*normrₓ*transpose(Sₓ)*rₓ
            ∂f∂eᵀ = 1/f*transpose(-rₓₓ×Sₓ .+ rₓ×Sₓₓ)*(rrₓₓₓ)
            ∂κ∂eᵀ = @. 1/g^2*(g*∂f∂eᵀ-f*∂g∂eᵀ)
            E*I*κ*∂κ∂eᵀ
        end
        QA = sum(w*qA(x) for (x,w) in zip(g5x,g5w))
        QI = sum(w*qI(x) for (x,w) in zip(g3x,g3w))
        QI .+ QA
    end
end

function make_V(ancs::ANC3DRURU)
    (;E,A,I) = ancs
    _x = find_x(ancs)
    g3ξ,g3w = FGQ.gausslegendre(3)
    g5ξ,g5w = FGQ.gausslegendre(5)
    g3x = _x.(g3ξ)
    g5x = _x.(g5ξ)
    function V(e)
        _rₓ  = make_rₓ(ancs,e)
        _rₓₓ = make_rₓₓ(ancs,e)   
        # _κ = make_κ(ancs,e)
        # _∂κ∂eᵀ = make_∂κ∂eᵀ(ancs,e)
        # _εₓₓ = make_εₓₓ(ancs,e)
        function qA(x)
            rₓ = _rₓ(x)
            εₓₓ = (transpose(rₓ)*rₓ - 1)/2
            εₓₓ^2
        end
        function qI(x) 
            rₓ = _rₓ(x)
            rₓₓ = _rₓₓ(x)
            normrₓ = norm(rₓ)
            rrₓₓₓ = rₓ × rₓₓ
            f = norm(rrₓₓₓ)
            g = normrₓ^3
            κ = f/g
            κ^2
        end
        VA = A.*sum(w*qA(x) for (x,w) in zip(g5x,g5w))
        VI = I.*sum(w*qI(x) for (x,w) in zip(g3x,g3w))
        E.*(VI .+ VA)./2
    end
end

function make_∂Q∂e(ancs::ANC3DRURU)
    Q = make_Q(ancs)
    function ∂Q∂e(e)
        ne = length(e)
        eT = eltype(e)
        out = zeros(eT,ne,ne)
        ∂Q∂e = ForwardDiff.jacobian!(out,Q,e)
    end
end

function build_M(ancs::ANC3DRURU)
    (;A,L,ρ) = ancs
    _S  = make_S(ancs)
    function integrand(x)
        S = _S(x)
        A.*ρ(x).*transpose(S)*S
    end
    integral, err = quadgk(integrand, 0, L)
    integral
end

function build_Sg(ancs::ANC3DRURU,mass)    
    (;A,L,ρ) = ancs
    _S  = make_S(ancs)
    function integrand(x)
        ρ(x).*_S(x)
    end
    integral, err = quadgk(integrand, 0, L)
    A.*integral./mass
end

function build_mass(ancs::ANC3DRURU)
    (;ρ) = ancs
    build_mass(ancs,ρ)
end

function build_mass(ancs::ANC3DRURU,ρ)
    (;A,L) = ancs
    function integrand(x)
        ρ(x)
    end
    integral, err = quadgk(integrand, 0, L)
    A*integral
end

function build_G(ancs::ANC3DRURU{T},g=9.81) where T
    f = SVector{3,T}(0,0,-g)
    mass = build_mass(ancs)
    Sg = build_Sg(ancs,mass)
    G = mass.*transpose(Sg)*f
    G
end

# CoordinateFunctions
"""
Encapsulated functions for ANC
$(TYPEDEF)
"""
struct CoordinateFunctions{ancsType,xT<:Function,ST<:Function}
    ancs::ancsType
    x::xT
    S::ST
end

function CoordinateFunctions(ancs,free_idx_mask,Φ_mask)
    _x = make_x(ancs)
    _S = make_S(ancs)
    CoordinateFunctions(ancs,_x,_S)
end


"""
返回未约束的Absolute Nodal Coordinates编号。
$(TYPEDSIGNATURES)
"""
function get_unconstrained_indices(ancs::ANC,pres_idx)
    deleteat!(collect(1:get_ncoords(ancs)),pres_idx)
end

"""
封装有函数的Absolute Nodal Coordinates类构造子。
$(TYPEDSIGNATURES)
"""
function CoordinateFunctions(ancs,uci)
    ξ = find_ξ(ancs)
    S = make_S(ancs)
    nq = get_ncoords(ancs)
    nuc = length(uci)
    CoordinateFunctions(ancs,ξ,S)
end

end
