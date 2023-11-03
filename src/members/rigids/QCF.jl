module QCF
using LinearAlgebra
using Rotations
using Quaternions
using StaticArrays
using ForwardDiff
using DocStringExtensions


import ..Rible: get_num_of_cstr, get_num_of_coords
import ..Rible: get_num_of_dof, get_num_of_local_dims
# import ..Rible: to_local_coords, to_transformation

import ..Rible: make_cstr_function
import ..Rible: make_cstr_jacobian
import ..Rible: make_cstr_hessians
import ..Rible: make_cstr_forces_jacobian
import ..Rible: cartesian_frame2coords

export get_num_of_cstr, get_num_of_coords
export get_num_of_dof, get_num_of_local_dims
# export to_local_coords, to_transformation

export make_cstr_function
export make_cstr_jacobian
export make_cstr_hessians
export make_cstr_forces_jacobian
export cartesian_frame2coords

vec(q::Quaternion) = SA[q.s, q.v1, q.v2, q.v3]

function Lmat(q::AbstractVector)
    q0,q1,q2,q3 = q
    SA[
        -q1  q0  q3 -q2;
        -q2 -q3  q0  q1;
        -q3  q2 -q1  q0;
    ]
end

function Lᵀmat(q::AbstractVector)
    q0,q1,q2,q3 = q
    SA[
        -q1 -q2 -q3;
         q0 -q3  q2;
         q3  q0 -q1;
        -q2  q1  q0;
    ]
end

function ∂Lᵀη∂q(η::AbstractVector)
    η1,η2,η3 = η
    SA[
         0  -η1 -η2 -η3;
        η1    0  η3 -η2;
        η2  -η3   0  η1;
        η3   η2 -η1   0;
    ]
end

function quatVel2localAngular(q::AbstractVector,q̇::AbstractVector)
    2Lmat(q)*q̇
end

function localAngular2quatVel(q::AbstractVector,Ω::AbstractVector)
    Lᵀmat(q)*Ω./2
end

function mass_matrix(J,q::AbstractVector)
    4Lᵀmat(q)*J*Lmat(q)
end

function get_Mγ(Jγ,γ,q::AbstractVector)
    4(Lᵀmat(q)*Jγ*Lmat(q) + γ*I)
end

function get_M⁻¹γ(J⁻¹γ,γ⁻¹,q::AbstractVector)
    (1//4)*(Lᵀmat(q)*J⁻¹γ*Lmat(q) + γ⁻¹*I)
end

function get_∂Mγq̇∂q(Jγ,q::AbstractVector,q̇::AbstractVector)
    L = Lmat(q)
    η = Jγ*L*q̇
    -4(transpose(L)*Jγ*Lmat(q̇) .+ ∂Lᵀη∂q(η))
end

function get_∂M⁻¹γp∂q(J⁻¹γ,q::AbstractVector,p::AbstractVector)
    L = Lmat(q)
    η = J⁻¹γ*L*p
    -1//4*(transpose(L)*J⁻¹γ*Lmat(p) .+ ∂Lᵀη∂q(η))
end

function get_∂Tγ∂qᵀ(Jγ,q::AbstractVector,q̇::AbstractVector)
    4Lᵀmat(q̇)*Jγ*Lmat(q̇)*q
end

function get_∂Tγ∂qᵀ∂q(Jγ,q̇::AbstractVector)
    4Lᵀmat(q̇)*Jγ*Lmat(q̇)
end

function Rmat(q::AbstractVector)
    q0,q1,q2,q3 = q
    v = SA[q1,q2,q3]
    q0² = q0^2
    q1² = q1^2
    q2² = q2^2
    q3² = q3^2
    Symmetric(SA[
        q0² + q1² - q2² - q3²     2q1*q2        2q1*q3;
           0        q0² - q1² + q2² - q3²       2q2*q3;
           0           0          q0² - q1² - q2² + q3²;
    ]) + 
    skew(2q0.*v)
end

function ∂Rη∂q(q::AbstractVector,η::AbstractVector)
    q0,q1,q2,q3 = q
    v = SA[q1,q2,q3]
    Q = q0*I + skew(v)
    Qη = Q*η    
    2hcat(
        Qη,
        transpose(v)*η*I - skew(Qη)
    )
end

function ∂Rᵀη∂q(q::AbstractVector,η::AbstractVector)
    q0,q1,q2,q3 = q
    v = SA[q1,q2,q3]
    Q = q0*I - skew(v)
    Qη = Q*η    
    2hcat(
        Qη,
        transpose(v)*η*I + skew(Qη)
    )
end

struct QC{T,JT}
    m::T
    m⁻¹::T
    γ::T
    Jγ::JT
    γ⁻¹::T
    J⁻¹γ::JT
end

function QC(m::T,J::AbstractMatrix{T};γ=maximum(diag(J))) where {T}
    m⁻¹ = inv(m)
    Jγ = Diagonal(J) - γ*I
    γ⁻¹ = inv(γ)
    J⁻¹γ = Diagonal(inv(J)) - γ⁻¹*I
    QC(m,m⁻¹,γ,Jγ,γ⁻¹,J⁻¹γ)
end

get_num_of_cstr(::QC) = 1
get_num_of_coords(::QC) = 7
get_num_of_dof(::QC) = 6
get_num_of_local_dims(::QC) = 3

function make_cstr_function()
    function inner_cstr_function(x::AbstractVector)
        q = @view x[4:7]
        (transpose(q)*q - 1)/2
    end
end

function make_cstr_jacobian()
    function inner_cstr_jacobian(x::AbstractVector)
        q = @view x[4:7]
        o = zero(eltype(q))
        SA[
            o o o q[1] q[2] q[3] q[4];
        ]
    end
end

"""
返回斜对称矩阵。
$(TYPEDSIGNATURES)
"""
function skew(w)
    w1,w2,w3 = w
    o = zero(w1)
    SA[
         o -w3 w2;
         w3 o -w1;
        -w2 w1 o
    ]
end

function to_transformation(c)
    ĉ = skew(c)
    function inner_C(x)
        q = @view x[4:7]
        T = eltype(q)
        I3 = SMatrix{3,3}(one(T)*I)
        hcat(
            I3,
            -2Rmat(q)*ĉ*Lmat(q)
        )
    end
end

function make_∂Cẋ∂x(c)
    ĉ = skew(c)
    function ∂Cẋ∂x(x,ẋ)
        q = @view x[4:7]
        q̇ = @view ẋ[4:7]
        R = Rmat(q)
        ĉL̇ = ĉ*Lmat(q̇)
        η = ĉL̇*q
        T = eltype(q)
        O3 = @SMatrix zeros(T,3,3)
        hcat(
            O3,
            2(R*ĉL̇ + ∂Rη∂q(q,η))
        )        
    end
end

function make_∂Cᵀf∂x(c)
    ĉᵀ = transpose(skew(c))
    function ∂Cᵀf∂x(x,f)
        q = @view x[4:7]
        ĉᵀRᵀf = ĉᵀ*transpose(Rmat(q))*f
        T = eltype(q)
        O37 = @SMatrix zeros(T,3,7)
        O43 = @SMatrix zeros(T,4,3)
        vcat(
            O37,
            hcat(
                O43,
                -2(Lᵀmat(q)*ĉᵀ*∂Rᵀη∂q(q,f) + ∂Lᵀη∂q(ĉᵀRᵀf))
            )
        )        
    end
end

function make_M(qcs::QC)
    (;m,γ,Jγ) = qcs
    function inner_M(x)
        q = @view x[4:7]
        T = eltype(x)
        O34 = @SMatrix zeros(T,3,4)
        hcat(
            vcat(
                SMatrix{3,3}(m*I),
                transpose(O34)
            ),
            vcat(
                O34,
                get_Mγ(Jγ,γ,q)
            )
        )
    end
end

function make_M⁻¹(qcs::QC)
    (;m⁻¹,γ⁻¹,J⁻¹γ) = qcs
    function inner_M(x)
        q = @view x[4:7]
        T = eltype(q)
        O34 = @SMatrix zeros(T,3,4)
        hcat(
            vcat(
                SMatrix{3,3}(m⁻¹*I),
                transpose(O34)
            ),
            vcat(
                O34,
                get_M⁻¹γ(J⁻¹γ,γ⁻¹,q)
            )
        )
    end
end

function make_∂Mẋ∂x(qcs::QC)
    (;Jγ) = qcs
    function inner_∂Mẋ∂x(x,ẋ)
        q = @view x[4:7]
        q̇ = @view ẋ[4:7]
        T = eltype(q)
        O37 = @SMatrix zeros(T,3,7)
        O43 = @SMatrix zeros(T,4,3)
        vcat(
            O37,
            hcat(
                O43,
                get_∂Mγq̇∂q(Jγ,q,q̇)
            )
        )
    end
end

function make_∂M⁻¹y∂x(qcs::QC)
    (;J⁻¹γ) = qcs
    function inner_∂M⁻¹y∂x(x,y)
        q = @view x[4:7]
        p = @view y[4:7]
        T = eltype(q)
        O37 = @SMatrix zeros(T,3,7)
        O43 = @SMatrix zeros(T,4,3)
        vcat(
            O37,
            hcat(
                O43,
                get_∂M⁻¹γp∂q(J⁻¹γ,q,p)
            )
        )
    end
end

function make_∂Cẋ∂x_forwarddiff(C,nc,nx)
    function ∂Cẋ∂x(x,ẋ)
        function Cẋ(_x)
            C(_x)*ẋ
        end
        out = zeros(eltype(ẋ),nc,nx)
        ForwardDiff.jacobian!(out,Cẋ,x)
    end
end

function make_∂T∂xᵀ(qcs::QC)
    (;Jγ) = qcs
    function inner_∂T∂xᵀ(x,ẋ)
        q = @view x[4:7]
        q̇ = @view ẋ[4:7]
        T = eltype(q)
        o3 = @SVector zeros(T,3)
        vcat(
            o3,
            get_∂Tγ∂qᵀ(Jγ,q,q̇)
        )            
    end
end

function make_∂T∂xᵀ∂x(qcs::QC)
    (;Jγ) = qcs
    function inner_∂T∂xᵀ∂x(ẋ)
        q̇ = @view ẋ[4:7]
        T = eltype(q)
        O37 = @SMatrix zeros(T,3,7)
        O43 = @SMatrix zeros(T,4,3)
        vcat(
            O37,
            hcat(
                O43,
                get_∂Tγ∂qᵀ∂q(Jγ,q̇)
            )
        )
    end
end

function cartesian_frame2coords(::QC,origin_position,R,origin_velocity,ω)
    Rmat = RotMatrix(R)
    q = QuatRotation(Rmat).q |> vec
    Ω = inv(Rmat)*ω
    q̇ = localAngular2quatVel(q,Ω)
    x = vcat(origin_position,q)
    ẋ = vcat(origin_velocity,q̇)
    x, ẋ
end

function find_rotation(x)
    q = @view x[4:7]
    RotMatrix(QuatRotation(q))
end

function find_local_angular_velocity(x,ẋ)
    q = @view x[4:7] 
    q̇ = @view ẋ[4:7]
    quatVel2localAngular(q,q̇)
end

function find_angular_velocity(x,ẋ)
    R = find_rotation(x)
    Ω = find_local_angular_velocity(x,ẋ)
    R*Ω
end


cstr_forces_jacobian(λ::AbstractVector) = cstr_forces_jacobian(first(λ))

function cstr_forces_jacobian(λ)
    o = zero(λ)    
    Diagonal(SA[o,o,o,λ,λ,λ,λ])
end

struct CoordinateFunctions{QCT,MT,M⁻¹T,∂Mẋ∂xT,∂M⁻¹y∂xT,∂T∂xᵀT,∂T∂xᵀ∂xT,cstr_forces_jacobianT,ΦT,ΦqT,cT}
    nmcs::QCT
    build_M::MT
    build_M⁻¹::M⁻¹T
    build_∂Mẋ∂x::∂Mẋ∂xT
    build_∂M⁻¹y∂x::∂M⁻¹y∂xT
    build_∂T∂xᵀ::∂T∂xᵀT
    build_∂T∂xᵀ∂x::∂T∂xᵀ∂xT
    cstr_forces_jacobian::cstr_forces_jacobianT
    Φ::ΦT
    Φq::ΦqT
    c::cT
end


function CoordinateFunctions(qcs)
    build_M = make_M(qcs)
    build_M⁻¹ = make_M⁻¹(qcs)
    build_∂Mẋ∂x = make_∂Mẋ∂x(qcs)
    build_∂M⁻¹y∂x = make_∂M⁻¹y∂x(qcs)
    build_∂T∂xᵀ = make_∂T∂xᵀ(qcs)
    build_∂T∂xᵀ∂x = make_∂T∂xᵀ∂x(qcs)    
    Φ = make_cstr_function()
    Φq = make_cstr_jacobian()
    c(x) = x
    CoordinateFunctions(
        qcs,
        build_M,
        build_M⁻¹,
        build_∂Mẋ∂x,
        build_∂M⁻¹y∂x,
        build_∂T∂xᵀ,
        build_∂T∂xᵀ∂x,
        cstr_forces_jacobian,
        Φ,
        Φq,
        c,
    )
end

end