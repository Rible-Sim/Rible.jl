q = (@SVector rand(4)) |> normalize
q̇ = (@SVector rand(4)) |> normalize
η = (@SVector rand(3)) |> normalize
Ω = (@SVector rand(3)) |> normalize

@btime RB.QCF.Lmat($q)
@btime RB.QCF.Lᵀmat($q)
@btime RB.QCF.∂Lᵀη∂q($η)
@btime RB.QCF.Rmat($q)
@btime RB.QCF.∂Rη∂q($q,$η)
@btime RB.QCF.∂Rᵀη∂q($q,$η)
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
