using Test
using Rible
import Rible as RB
using RibleQCF
import RibleQCF as QCF
using Random
using StaticArrays
using ForwardDiff
using LinearAlgebra

@testset "QCF" begin
    q = (@SVector rand(4)) |> normalize
    q̇ = (@SVector rand(4)) |> normalize
    η = (@SVector rand(3)) |> normalize
    Ω = (@SVector rand(3)) |> normalize

    @allocations QCF.Lmat(q)
    @allocations QCF.Lᵀmat(q)
    @allocations QCF.∂Lᵀη∂q(η)
    @allocations QCF.Rmat(q)
    @allocations QCF.∂Rη∂q(q,η)
    @allocations QCF.∂Rᵀf∂q(q,η)
    @allocations QCF.quatVel2localAngular(q, q̇)
    @allocations QCF.localAngular2quatVel(q, Ω)

    x = (@MVector rand(7)) 
    normalize!(@view x[4:7])
    m = rand()
    J = rand(3,3)
    qcs = QCF.QC(m,J)
    @allocations QCF.find_rotation(qcs,x)
    (;γ,Jγ,γ⁻¹,J⁻¹γ) = qcs
    Mγ = QCF.get_Mγ(Jγ,γ,q)
    @allocations QCF.get_Mγ(Jγ,γ,q)
    p = Mγ*q̇
    @allocations QCF.get_M⁻¹γ(J⁻¹γ,γ⁻¹,q)
    @allocations QCF.get_∂Mγq̇∂q(Jγ,q,q̇)
    @allocations QCF.get_∂M⁻¹γp∂q(J⁻¹γ,q,p)
    @allocations QCF.get_∂Tγ∂qᵀ(Jγ,q,q̇)
    @allocations QCF.get_∂Tγ∂qᵀ∂q(Jγ,q̇)


    η = @SVector rand(3)
    q = normalize(@SVector rand(4))
    ω = @SVector rand(3)
    q̇ = QCF.Lᵀmat(q)*ω
    QCF.∂Rη∂q(q,η)*q̇
    -2QCF.Rmat(q)*RB.skew(η)*QCF.Lmat(q)*q̇

    ro = @SVector rand(3)
    x = vcat(ro,q)
    qcs = QCF.QC(1.0,rand(3,3))

    @test has_constant_mass_matrix(qcs) === Val(false)

    c = @SVector rand(3)
    ṙo = @SVector rand(3)
    ẋ = vcat(ṙo,q̇)
    QCF.to_position_jacobian(qcs,x,c)

    ∂Cẋ∂x_ref = ForwardDiff.jacobian((x) -> QCF.to_position_jacobian(qcs,x,c)*ẋ,x)
    to_velocity_jacobian = QCF.to_velocity_jacobian(qcs, x,ẋ,c)
    @test to_velocity_jacobian ≈ ∂Cẋ∂x_ref

    ∂²Rη∂qᵀ∂q_ref = ForwardDiff.jacobian((q) -> QCF.∂Rη∂q(q,c),q)
    QCF.∂²Rη∂qᵀ∂q(c)

    reshape(∂²Rη∂qᵀ∂q_ref,3,4,4)
    H = QCF.∂²Rη∂qᵀ∂q(η)
end
