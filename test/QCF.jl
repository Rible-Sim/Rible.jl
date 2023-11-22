q = (@SVector rand(4)) |> normalize
q̇ = (@SVector rand(4)) |> normalize
η = (@SVector rand(3)) |> normalize
Ω = (@SVector rand(3)) |> normalize

@btime RB.QCF.Lmat($q)
@btime RB.QCF.Lᵀmat($q)
@btime RB.QCF.∂Lᵀη∂q($η)
@btime RB.QCF.Rmat($q)
@btime RB.QCF.∂Rη∂q($q,$η)
@btime RB.QCF.∂Rᵀf∂q($q,$η)
@btime RB.QCF.quatVel2localAngular($q, $q̇)
@btime RB.QCF.localAngular2quatVel($q, $Ω)

x = (@MVector rand(7)) 
normalize!(@view x[4:7])
m = rand()
J = rand(3,3)
qcs = RB.QCF.QC(m,J)
@btime RB.QCF.find_rotation($qcs,$x)
(;γ,Jγ,γ⁻¹,J⁻¹γ) = qcs
Mγ = RB.QCF.get_Mγ(Jγ,γ,q)
@btime RB.QCF.get_Mγ($Jγ,$γ,$q)
p = Mγ*q̇
@btime RB.QCF.get_M⁻¹γ($J⁻¹γ,$γ⁻¹,$q)
@btime RB.QCF.get_∂Mγq̇∂q($Jγ,$q,$q̇)
@btime RB.QCF.get_∂M⁻¹γp∂q($J⁻¹γ,$q,$p)
@btime RB.QCF.get_∂Tγ∂qᵀ($Jγ,$q,$q̇)
@btime RB.QCF.get_∂Tγ∂qᵀ∂q($Jγ,$q̇)



η = @SVector rand(3)
q = normalize(@SVector rand(4))
ω = @SVector rand(3)
q̇ = RB.QCF.Lᵀmat(q)*ω
RB.QCF.∂Rη∂q(q,η)*q̇
-2RB.QCF.Rmat(q)*RB.skew(η)*RB.QCF.Lmat(q)*q̇

ro = @SVector rand(3)
x = vcat(ro,q)
qcs = RB.QCF.QC(1.0,rand(3,3))
c = @SVector rand(3)
ṙo = @SVector rand(3)
ẋ = vcat(ṙo,q̇)
RB.QCF.to_transformation(qcs,x,c)

∂Cẋ∂x_ref = ForwardDiff.jacobian((x) -> RB.QCF.to_transformation(qcs,x,c)*ẋ,x)
∂Cẋ∂x = RB.QCF.make_∂Cẋ∂x(c)(x,ẋ)
∂Cẋ∂x - ∂Cẋ∂x_ref

∂²Rη∂qᵀ∂q_ref = ForwardDiff.jacobian((q) -> RB.QCF.∂Rη∂q(q,c),q)
RB.QCF.∂²Rη∂qᵀ∂q(c)

reshape(∂²Rη∂qᵀ∂q_ref,3,4,4)
H = RB.QCF.∂²Rη∂qᵀ∂q(η)
