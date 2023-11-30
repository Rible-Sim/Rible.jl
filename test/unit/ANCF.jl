using Rible

@testset "ANCF" begin
    ancs = ANCF.ANC3DRURU(8.96e3;E=110e6,L,radius=0.01)
    S = ANCF.make_S(ancs)(0.5)
    Sₓ = ANCF.make_Sₓ(ancs)(0.5)
    Sₓₓ = ANCF.make_Sₓₓ(ancs)(0.5)
    ri = rand(3)
    riₓ = rand(3)
    rj = rand(3)
    rjₓ = rand(3)
    e = [ri;riₓ;rj;rjₓ]
    r = ANCF.make_r(ancs,e)(0.5)
    rₓ = ANCF.make_rₓ(ancs,e)(0.5)
    rₓₓ = ANCF.make_rₓₓ(ancs,e)(0.5)
    κ = ANCF.make_κ(ancs,e)(0.5)
    _κ = ANCF.make_κ(ancs,e)
    function ∂κ∂eᵀ_forwarddiff!(out,ancs,x,e)
        function κ(e)
            ANCF.make_κ(ancs,e)(x)
        end
        ForwardDiff.gradient!(out,κ,e)
    end
    ne = length(e)
    eT = eltype(e)
    out = zeros(eT,ne)
    ∂κ∂eᵀ_forwarddiff!(out,ancs,0.5L,e)
    @test ANCF.make_∂κ∂eᵀ(ancs,e)(0.5) ≈ ∂κ∂eᵀ_forwarddiff!(out,ancs,0.5,e)

    ∂κ∂eᵀ_forwarddiff!(out,ancs,0.5,e)

    ANCF.make_S(ancs)(0.5)
    ANCF.make_Sₓ(ancs)(0.5)
    ANCF.make_Sₓₓ(ancs)(0.5)
    ANCF.make_r(ancs,e)(0.5)
    ANCF.make_rₓ(ancs,e)(0.5)
    ANCF.make_rₓₓ(ancs,e)(0.5)
    ANCF.make_κ(ancs,e)(0.5)
    ANCF.make_∂κ∂eᵀ(ancs,e)(0.5)
    Q = ANCF.make_Q(ancs)

    ANCF.build_M(ancs)

    M = ANCF.build_M(ancs)
    G = ANCF.build_G(ancs) 

    # ne = length(e)
    # eT = eltype(e)
    # out = zeros(eT,ne)
    # ∂Q∂e_forwarddiff!(out,ancs,0.5L,e)
    # ForwardDiff.jacobian!(out,Q,e)
    # ∂Q∂e |> issymmetric
end
