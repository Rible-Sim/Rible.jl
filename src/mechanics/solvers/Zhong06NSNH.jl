
function Momentum_k(qᵏ⁻¹,pᵏ⁻¹,qᵏ,λᵏ,μᵏ,M,A,B,h)
    pᵏ = -pᵏ⁻¹ .+ 2/h.*M*(qᵏ.-qᵏ⁻¹) .+
        1/(2h).*(transpose(A(qᵏ))-transpose(A(qᵏ⁻¹)))*λᵏ .+
        1/(2h).*(transpose(B(qᵏ))-transpose(B(qᵏ⁻¹)))*μᵏ
end

function stepk_maker(nq,nλ,nμ,qₛ₋₁,q̇ₛ₋₁,pₛ₋₁,tₛ₋₁,dynfuncs,invM,h)
    M,Φ,A,Ψ,B,F!,jacobians,contact_funcs = dynfuncs
    Jac_F!,Ψq,∂Aᵀλ∂q,∂Bᵀμ∂q = jacobians
    # E,𝐠,𝐠𝐪,∂𝐠𝐪ᵀΛ∂q,∂𝐠𝐪q̇∂q = contact_funcs
    n1 =  nq
    n2 = nq+nλ
    n3 = nq+nλ+nμ
    function stepk!(R,J,x,q̇ₛ)
        qₛ = @view x[   1:n1]
        λₛ = @view x[n1+1:n2]
        μₛ = @view x[n2+1:n3]

        q = (qₛ.+qₛ₋₁)./2
        q̇ = (qₛ.-qₛ₋₁)./h
        t = tₛ₋₁+h/2

        F⁺ = zeros(eltype(qₛ),nq)
        F!(F⁺,q,q̇,t)
        ∂F∂q = zeros(eltype(qₛ),nq,nq)
        ∂F∂q̇ = zeros(eltype(qₛ),nq,nq)
        Jac_F!(∂F∂q,∂F∂q̇,q,q̇,t)

        # q̇ₛ = invM*Momentum_k(qₛ₋₁,pₛ₋₁,qₛ,λₛ,μₛ,M,A,B,h)

        Aₛ₋₁ = A(qₛ₋₁)
        Bₛ₋₁ = B(qₛ₋₁)
        Aₛ = A(qₛ)
        Bₛ = B(qₛ)
        ∂q̇ₛ∂qₛ = 2/h*I + 1/(2h).*invM*(∂Aᵀλ∂q(qₛ,λₛ) + ∂Bᵀμ∂q(qₛ,μₛ))
        ∂q̇ₛ∂λₛ = invM*transpose(Aₛ-Aₛ₋₁)/(2h)
        ∂q̇ₛ∂μₛ = invM*transpose(Bₛ-Bₛ₋₁)/(2h)

        R[   1:n1] .= -h.*pₛ₋₁ .+ M*(qₛ.-qₛ₋₁) .-
                        1/2 .*transpose(Aₛ₋₁)*λₛ .-
                        1/2 .*transpose(Bₛ₋₁)*μₛ .-
                        (h^2)/2 .*F⁺
        R[n1+1:n2] .= Φ(qₛ)
        R[n2+1:n3] .= Ψ(qₛ,q̇ₛ)

        J .= 0.0
        J[   1:n1,   1:n1] .=  M .-h^2/2 .*(1/2 .*∂F∂q .+ 1/h.*∂F∂q̇)
        J[   1:n1,n1+1:n2] .= -1/2 .*transpose(Aₛ₋₁)
        J[   1:n1,n2+1:n3] .= -1/2 .*transpose(Bₛ₋₁)

        J[n1+1:n2,   1:n1] .=  Aₛ

        J[n2+1:n3,   1:n1] .=  Ψq(qₛ,q̇ₛ)+Bₛ*∂q̇ₛ∂qₛ
        J[n2+1:n3,n1+1:n2] .=  Bₛ*∂q̇ₛ∂λₛ
        J[n2+1:n3,n2+1:n3] .=  Bₛ*∂q̇ₛ∂μₛ
    end

    stepk!
end

function ns_stepk_maker(nq,nλ,nμ,nu,qᵏ⁻¹,q̇ᵏ⁻¹,pᵏ⁻¹,tᵏ⁻¹,Aset,dynfuncs,invM,h)
    M,Φ,A,Ψ,B,F!,jacobians,contact_funcs = dynfuncs
    Jac_F!,Ψq,∂Aᵀλ∂q,∂Bᵀμ∂q = jacobians
    E,𝐠,𝐠𝐪,∂𝐠𝐪ᵀΛ∂q,∂𝐠𝐪q̇∂q = contact_funcs
    r = 1.0
    n1 =  nq
    n2 = 2nq
    n3 = 2nq+nλ
    n4 = 2nq+nλ+nμ
    function inner_ns_stepk!(R,J,x,q̇ᵏ)
        q̃ᵏ = @view x[   1:n1]
        qᵏ = @view x[n1+1:n2]
        λᵏ = @view x[n2+1:n3]
        μᵏ = @view x[n3+1:n4]
        Λᵏ = @view x[n4+1:n4+nu]

        q = (qᵏ.+qᵏ⁻¹)./2
        q̇ = (qᵏ.-qᵏ⁻¹)./h
        t = tᵏ⁻¹+h/2

        F⁺ = zeros(eltype(qᵏ),nq)
        F!(F⁺,q,q̇,t)
        ∂F∂q = zeros(eltype(qᵏ),nq,nq)
        ∂F∂q̇ = zeros(eltype(qᵏ),nq,nq)
        Jac_F!(∂F∂q,∂F∂q̇,q,q̇,t)

        Aᵏ⁻¹ = A(qᵏ⁻¹)
        Bᵏ⁻¹ = B(qᵏ⁻¹)
        Aᵏ = A(qᵏ)
        Bᵏ = B(qᵏ)
        ∂q̇ᵏ∂qᵏ = 1/(2h).*invM*(∂Aᵀλ∂q(qᵏ,λᵏ) .+ ∂Bᵀμ∂q(qᵏ,μᵏ)) + 2/h*I
        ∂q̇ᵏ∂λᵏ = invM*transpose(Aᵏ-Aᵏ⁻¹)/(2h)
        ∂q̇ᵏ∂μᵏ = invM*transpose(Bᵏ-Bᵏ⁻¹)/(2h)

        gqᵏ⁻¹ = 𝐠𝐪(qᵏ⁻¹)
        gqᵏ = 𝐠𝐪(qᵏ)
        Bset = Aset .& (Λᵏ .- r.*(gqᵏ*q̇ᵏ .+ E*gqᵏ⁻¹*q̇ᵏ⁻¹) .≥ 0)
        B̄set = .!Bset
        # @show q̇ᵏ[1], qᵏ[1]
        # if count(Aset)!=0
        #     @show gqᵏ*q̇ᵏ .+ E*gqᵏ⁻¹*q̇ᵏ⁻¹
        #     @show Λᵏ[1]
        #     @show Aset, Bset, B̄set
        # end
        nB = count(Bset)
        nB̄ = count(B̄set)
        Bindex = findall(Bset)
        B̄index = findall(B̄set)
        # if nB != 0
        #     @show Bindex,B̄index
        # end
        gqᵏB = gqᵏ[Bindex,:]
        gqᵏ⁻¹B = gqᵏ⁻¹[Bindex,:]
        n5 = 2nq+nλ+nμ+nB
        n6 = 2nq+nλ+nμ+nu

        R .= 0.0
        R[   1:n1] .= -h.*pᵏ⁻¹ .+ M*(q̃ᵏ.-qᵏ⁻¹) .-
                        1/2 .*transpose(Aᵏ⁻¹)*λᵏ .-
                        1/2 .*transpose(Bᵏ⁻¹)*μᵏ .-
                        (h^2)/2 .*F⁺
        R[n1+1:n2] .= (2/h).*M*(qᵏ - q̃ᵏ) - transpose(gqᵏ)*Λᵏ
        R[n2+1:n3] .= Φ(qᵏ)
        R[n3+1:n4] .= Ψ(qᵏ,q̇ᵏ)
        R[n4+1:n5] .= gqᵏB*q̇ᵏ .+ E[Bindex,Bindex]*gqᵏ⁻¹B*q̇ᵏ⁻¹
        R[n5+1:n6] .= Λᵏ[B̄index]
        R1 = -h.*pᵏ⁻¹ .+ M*(qᵏ.-qᵏ⁻¹) .-
                (h/2).*transpose(gqᵏ)*Λᵏ .-
                1/2 .*transpose(Aᵏ⁻¹)*λᵏ .-
                1/2 .*transpose(Bᵏ⁻¹)*μᵏ .-
                (h^2)/2 .*F⁺
        @show norm(R1,Inf),norm(R[   1:n1],Inf), norm(R[n1+1:n2],Inf)
        J .= 0.0
        J[   1:n1,   1:n1] .=  M
        J[   1:n1,n1+1:n2] .= -h^2/2 .*(1/2 .*∂F∂q .+ 1/h.*∂F∂q̇)
        J[   1:n1,n2+1:n3] .= -1/2 .*transpose(Aᵏ⁻¹)
        J[   1:n1,n3+1:n4] .= -1/2 .*transpose(Bᵏ⁻¹)

        J[n1+1:n2,   1:n1] .= -(2/h).*M
        J[n1+1:n2,n1+1:n2] .=  (2/h).*M .- ∂𝐠𝐪ᵀΛ∂q(qᵏ,Λᵏ)
        J[n1+1:n2,n4+1:n5] .= -transpose(gqᵏB)

        J[n2+1:n3,n1+1:n2] .=  Aᵏ

        J[n3+1:n4,n1+1:n2] .=  Ψq(qᵏ,q̇ᵏ)+Bᵏ*∂q̇ᵏ∂qᵏ
        J[n3+1:n4,n2+1:n3] .=  Bᵏ*∂q̇ᵏ∂λᵏ
        J[n3+1:n4,n3+1:n4] .=  Bᵏ*∂q̇ᵏ∂μᵏ

        J[n4+1:n5,n1+1:n2] .=  ∂𝐠𝐪q̇∂q(qᵏ,q̇ᵏ)[Bindex,:]+gqᵏB*∂q̇ᵏ∂qᵏ
        J[n4+1:n5,n2+1:n3] .=  gqᵏB*∂q̇ᵏ∂λᵏ
        J[n4+1:n5,n3+1:n4] .=  gqᵏB*∂q̇ᵏ∂μᵏ

        J[n5+1:n6,n5+1:n6] .=  Matrix(1I,nB̄,nB̄)

        nB, Bindex, B̄index
    end
end

function sns_stepk_maker(nq,nλ,nμ,nu,qᵏ⁻¹,q̇ᵏ⁻¹,pᵏ⁻¹,tᵏ⁻¹,Aset,dynfuncs,invM,h)
    M,Φ,A,Ψ,B,F!,jacobians,contact_funcs = dynfuncs
    Jac_F!,Ψq,∂Aᵀλ∂q,∂Bᵀμ∂q = jacobians
    E,𝐠,𝐠𝐪,∂𝐠𝐪ᵀΛ∂q,∂𝐠𝐪q̇∂q = contact_funcs
    r = 1.0
    n1 = nq
    n2 = nq+nλ
    n3 = nq+nλ+nμ
    function inner_ns_stepk!(R,J,x,q̇ᵏ)
        qᵏ = @view x[   1:n1]
        λᵏ = @view x[n1+1:n2]
        μᵏ = @view x[n2+1:n3]
        Λᵏ = @view x[n3+1:n3+nu]

        q = (qᵏ.+qᵏ⁻¹)./2
        q̇ = (qᵏ.-qᵏ⁻¹)./h
        t = tᵏ⁻¹+h/2

        F⁺ = zeros(eltype(qᵏ),nq)
        F!(F⁺,q,q̇,t)
        ∂F∂q = zeros(eltype(qᵏ),nq,nq)
        ∂F∂q̇ = zeros(eltype(qᵏ),nq,nq)
        Jac_F!(∂F∂q,∂F∂q̇,q,q̇,t)

        Aᵏ⁻¹ = A(qᵏ⁻¹)
        Bᵏ⁻¹ = B(qᵏ⁻¹)
        Aᵏ = A(qᵏ)
        Bᵏ = B(qᵏ)
        ∂q̇ᵏ∂qᵏ = 1/h.*invM*(∂Aᵀλ∂q(qᵏ,λᵏ) .+ ∂Bᵀμ∂q(qᵏ,μᵏ)) + 2/h*I
        ∂q̇ᵏ∂λᵏ = invM*transpose(Aᵏ-Aᵏ⁻¹)/h
        ∂q̇ᵏ∂μᵏ = invM*transpose(Bᵏ-Bᵏ⁻¹)/h

        gqᵏ⁻¹ = 𝐠𝐪(qᵏ⁻¹)
        gqᵏ = 𝐠𝐪(qᵏ)
        Bset = Aset .& (Λᵏ .- r.*(gqᵏ*q̇ᵏ .+ E*gqᵏ⁻¹*q̇ᵏ⁻¹) .≥ 0)
        B̄set = .!Bset
        # @show q̇ᵏ[1], qᵏ[1]
        # if count(Aset)!=0
        #     @show gqᵏ*q̇ᵏ .+ E*gqᵏ⁻¹*q̇ᵏ⁻¹
        #     @show Λᵏ[1]
        #     @show Aset, Bset, B̄set
        # end
        nB = count(Bset)
        nB̄ = count(B̄set)
        Bindex = findall(Bset)
        B̄index = findall(B̄set)
        # if nB != 0
        #     @show Bindex,B̄index
        # end
        gqᵏB = gqᵏ[Bindex,:]
        gqᵏ⁻¹B = gqᵏ⁻¹[Bindex,:]
        n4 = nq+nλ+nμ+nB
        n5 = nq+nλ+nμ+nu

        R .= 0.0
        R[   1:n1] .= -h.*pᵏ⁻¹ .+ M*(qᵏ.-qᵏ⁻¹) .-
                        transpose(gqᵏ⁻¹)*Λᵏ .-
                        transpose(Aᵏ⁻¹)*λᵏ .-
                        transpose(Bᵏ⁻¹)*μᵏ .-
                        (h^2)/2 .*F⁺
        R[n1+1:n2] .= Φ(qᵏ)
        R[n2+1:n3] .= Ψ(qᵏ,q̇ᵏ)
        R[n3+1:n4] .= gqᵏB*q̇ᵏ .+ E[Bindex,Bindex]*gqᵏ⁻¹B*q̇ᵏ⁻¹
        R[n4+1:n5] .= Λᵏ[B̄index]
        J .= 0.0
        J[   1:n1,   1:n1] .=  M .- h^2/2 .*(1/2 .*∂F∂q .+ 1/h.*∂F∂q̇) #.- ∂𝐠𝐪ᵀΛ∂q(qᵏ,Λᵏ)
        J[   1:n1,n1+1:n2] .= -transpose(Aᵏ⁻¹)
        J[   1:n1,n2+1:n3] .= -transpose(Bᵏ⁻¹)
        J[   1:n1,n3+1:n4] .= -transpose(gqᵏ⁻¹B)

        J[n1+1:n2,   1:n1] .=  Aᵏ

        J[n2+1:n3,   1:n1] .=  Ψq(qᵏ,q̇ᵏ)+Bᵏ*∂q̇ᵏ∂qᵏ
        J[n2+1:n3,n1+1:n2] .=  Bᵏ*∂q̇ᵏ∂λᵏ
        J[n2+1:n3,n2+1:n3] .=  Bᵏ*∂q̇ᵏ∂μᵏ

        J[n3+1:n4,   1:n1] .=  ∂𝐠𝐪q̇∂q(qᵏ,q̇ᵏ)[Bindex,:]+gqᵏB*∂q̇ᵏ∂qᵏ
        J[n3+1:n4,n1+1:n2] .=  gqᵏB*∂q̇ᵏ∂λᵏ
        J[n3+1:n4,n2+1:n3] .=  gqᵏB*∂q̇ᵏ∂μᵏ

        J[n4+1:n5,n4+1:n5] .=  Matrix(1I,nB̄,nB̄)

        nB, Bindex, B̄index
    end
end

function find_step_length(y,Λ,Δy,ΔΛ,η)
    αPindex = findall((x)->x<0,Δy)
    αDindex = findall((x)->x<0,ΔΛ)
    if isempty(αPindex)
        αP = one(eltype(y))
    else
        αP = min(one(eltype(y)),η*minimum(-y[αPindex]./Δy[αPindex]))
    end
    if isempty(αDindex)
        αD = one(eltype(Λ))
    else
        αD = min(one(eltype(Λ)),η*minimum(-Λ[αDindex]./ΔΛ[αDindex]))
    end
    αP, αD
end

function ip_ns_stepk_maker(nq,nλ,nμ,nu,qᵏ⁻¹,q̇ᵏ⁻¹,pᵏ⁻¹,tᵏ⁻¹,Aset,dynfuncs,invM,h)
    M,Φ,A,Ψ,B,F!,jacobians,contact_funcs = dynfuncs
    Jac_F!,Ψq,∂Aᵀλ∂q,∂Bᵀμ∂q = jacobians
    E,𝐠,𝐠𝐪,∂𝐠𝐪ᵀΛ∂q,∂𝐠𝐪q̇∂q = contact_funcs
    r = 1.0
    n1 = nq
    n2 = n1 + nq
    n3 = n2 + nλ
    n4 = n3 + nμ
    n5 = n4 + nu
    n6 = n5 + nu
    function inner_ip_ns_stepk!(R,J,x,q̇ᵏ,μ)
        q̃ᵏ = @view x[   1:n1]
        qᵏ = @view x[n1+1:n2]
        λᵏ = @view x[n2+1:n3]
        μᵏ = @view x[n3+1:n4]
        Λᵏ = @view x[n4+1:n5]
        yᵏ = @view x[n5+1:n6]

        q = (qᵏ.+qᵏ⁻¹)./2
        q̇ = (qᵏ.-qᵏ⁻¹)./h
        t = tᵏ⁻¹+h/2

        F⁺ = zeros(eltype(qᵏ),nq)
        F!(F⁺,q,q̇,t)
        ∂F∂q = zeros(eltype(qᵏ),nq,nq)
        ∂F∂q̇ = zeros(eltype(qᵏ),nq,nq)
        Jac_F!(∂F∂q,∂F∂q̇,q,q̇,t)

        Aᵏ⁻¹ = A(qᵏ⁻¹)
        Bᵏ⁻¹ = B(qᵏ⁻¹)
        Aᵏ = A(qᵏ)
        Bᵏ = B(qᵏ)
        ∂q̇ᵏ∂qᵏ = 1/(2h).*invM*(∂Aᵀλ∂q(qᵏ,λᵏ) .+ ∂Bᵀμ∂q(qᵏ,μᵏ)) + 2/h*I
        ∂q̇ᵏ∂λᵏ = invM*transpose(Aᵏ-Aᵏ⁻¹)/(2h)
        ∂q̇ᵏ∂μᵏ = invM*transpose(Bᵏ-Bᵏ⁻¹)/(2h)

        gqᵏ⁻¹ = 𝐠𝐪(qᵏ⁻¹)
        gqᵏ = 𝐠𝐪(qᵏ)
        Aindex = findall(Aset)
        gqᵏA = gqᵏ[Aindex,:]
        gqᵏ⁻¹A = gqᵏ⁻¹[Aindex,:]

        R .= 0.0
        R[   1:n1] .= -h.*pᵏ⁻¹ .+ M*(q̃ᵏ.-qᵏ⁻¹) .-
                        1/2 .*transpose(Aᵏ⁻¹)*λᵏ .-
                        1/2 .*transpose(Bᵏ⁻¹)*μᵏ .-
                        (h^2)/2 .*F⁺
        R[n1+1:n2] .= (2/h).*M*(qᵏ - q̃ᵏ) - transpose(gqᵏA)*Λᵏ
        R[n2+1:n3] .= Φ(qᵏ)
        R[n3+1:n4] .= Ψ(qᵏ,q̇ᵏ)
        R[n4+1:n5] .= gqᵏA*q̇ᵏ .+ E[Aindex,Aindex]*gqᵏ⁻¹A*q̇ᵏ⁻¹ .- yᵏ
        R[n5+1:n6] .= Λᵏ.*yᵏ

        J .= 0.0
        J[   1:n1,   1:n1] .=  M
        J[   1:n1,n1+1:n2] .= -h^2/2 .*(1/2 .*∂F∂q .+ 1/h.*∂F∂q̇)
        J[   1:n1,n2+1:n3] .= -1/2 .*transpose(Aᵏ⁻¹)
        J[   1:n1,n3+1:n4] .= -1/2 .*transpose(Bᵏ⁻¹)

        J[n1+1:n2,   1:n1] .= -(2/h).*M
        J[n1+1:n2,n1+1:n2] .=  (2/h).*M .- ∂𝐠𝐪ᵀΛ∂q(qᵏ,Λᵏ)
        J[n1+1:n2,n4+1:n5] .= -transpose(gqᵏA)

        J[n2+1:n3,n1+1:n2] .=  Aᵏ

        J[n3+1:n4,n1+1:n2] .=  Ψq(qᵏ,q̇ᵏ)+Bᵏ*∂q̇ᵏ∂qᵏ
        J[n3+1:n4,n2+1:n3] .=  Bᵏ*∂q̇ᵏ∂λᵏ
        J[n3+1:n4,n3+1:n4] .=  Bᵏ*∂q̇ᵏ∂μᵏ

        J[n4+1:n5,n1+1:n2] .=  ∂𝐠𝐪q̇∂q(qᵏ,q̇ᵏ)[Aindex,:]+gqᵏA*∂q̇ᵏ∂qᵏ
        J[n4+1:n5,n2+1:n3] .=  gqᵏA*∂q̇ᵏ∂λᵏ
        J[n4+1:n5,n3+1:n4] .=  gqᵏA*∂q̇ᵏ∂μᵏ
        J[n4+1:n5,n5+1:n6] .= -Matrix(1I,nu,nu)

        J[n5+1:n6,n4+1:n5] .=  Diagonal(yᵏ)
        J[n5+1:n6,n5+1:n6] .=  Diagonal(Λᵏ)

        η = 1.0
        Δxp = J\-R
        ΔΛp = @view Δxp[n4+1:n5]
        Δyp = @view Δxp[n5+1:n6]
        αP, αD = find_step_length(yᵏ,Λᵏ,Δyp,ΔΛp,η)
        yp = yᵏ .+ αP.*Δyp
        Λp = Λᵏ .+ αD.*ΔΛp
        μp = transpose(yp)*Λp/nu
        σ = (μp/μ)^3
        τ = σ*μ
        @show τ,σ,μ
        R[n5+1:n6] .+ Δyp .* ΔΛp .- σ.*μ
        Δxc = J\-R
        # η = exp(-0.1μ) + 0.9
        ΔΛc = @view Δxc[n4+1:n5]
        Δyc = @view Δxc[n5+1:n6]
        αP, αD = find_step_length(yᵏ,Λᵏ,Δyc,ΔΛc,η)
        q̃ᵏ .+= αP.*Δxc[   1:n1]
        qᵏ .+= αP.*Δxc[n1+1:n2]
        λᵏ .+= αD.*Δxc[n2+1:n3]
        μᵏ .+= αD.*Δxc[n3+1:n4]
        Λᵏ .+= αD.*ΔΛc
        yᵏ .+= αP.*Δyc
        μ = transpose(yᵏ)*Λᵏ/nu
    end
end


function nhsolve(prob,nq,nλ,nμ,nu,q0,q̇0;tspan,dt=0.01,ftol=1e-14,verbose=false,iterations=50,
                progress=true,exception=true)
    @unpack bot,dynfuncs = prob
    M,Φ,A,Ψ,B,F!,jacobians,contact_funcs = dynfuncs
    Jac_F!,Ψq,∂Aᵀλ∂q,∂Bᵀμ∂q = jacobians
    E,𝐠,𝐠𝐪,∂𝐠𝐪ᵀΛ∂q,∂𝐠𝐪q̇∂q = contact_funcs
    totaltime = tspan[end] - tspan[begin]
    totalstep = ceil(Int,totaltime/dt)
    ts = [tspan[begin]+(i-1)*dt for i in 1:totalstep+1]
    q̃s = [copy(q0) for i in 1:totalstep+1]
    qs = [copy(q0) for i in 1:totalstep+1]
    q̇s = [copy(q̇0) for i in 1:totalstep+1]
    ps = [M*copy(q̇0) for i in 1:totalstep+1]
    λs = [zeros(eltype(q0),nλ) for i in 1:totalstep+1]
    friction_coefficients = [zeros(eltype(q0),nμ) for i in 1:totalstep+1]
    Λs = [zeros(eltype(q0),nu) for i in 1:totalstep+1]
    invM = inv(M)
    # F⁺ = zero(q0)
    step = 0
    nx = nq + nq + nλ + nμ + nu
    initial_x = zeros(eltype(q0),nx)
    initial_R = zeros(eltype(q0),nx)
    initial_J = zeros(eltype(q0),nx,nx)
    smooth_nx = nq + nλ + nμ
    smooth_x = zeros(eltype(q0),smooth_nx)
    smooth_R = zeros(eltype(q0),smooth_nx)
    smooth_J = zeros(eltype(q0),smooth_nx,smooth_nx)
    mr = norm(M,Inf)
    scaling = 1

    iteration = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = 1:totalstep
        #---------Step k Control-----------
        # control!(intor,cache)
        #---------Step k Control-----------
        q̃ᵏ⁻¹ = q̃s[timestep]
        qᵏ⁻¹ = qs[timestep]
        q̇ᵏ⁻¹ = q̇s[timestep]
        pᵏ⁻¹ = ps[timestep]
        λᵏ⁻¹ = λs[timestep]
        μᵏ⁻¹ = friction_coefficients[timestep]
        Λᵏ⁻¹ = Λs[timestep]
        tᵏ⁻¹ = ts[timestep]
        q̃ᵏ   = q̃s[timestep+1]
        qᵏ   = qs[timestep+1]
        q̇ᵏ   = q̇s[timestep+1]
        pᵏ   = ps[timestep+1]
        λᵏ   = λs[timestep+1]
        μᵏ   = friction_coefficients[timestep+1]
        Λᵏ   = Λs[timestep+1]
        q̇ᵏ .= q̇ᵏ⁻¹

        qˣ = qᵏ⁻¹ .+ dt./2 .*q̇ᵏ⁻¹
        Aset = 𝐠(qˣ) .< 0
        Āset = .!Aset
        isconverged = false
        res = typemax(eltype(qᵏ))
        k_break = 0
        if count(Aset) == 0
            smooth_x[      1:nq]          .= qᵏ⁻¹
            smooth_x[   nq+1:nq+nλ]       .= 0.0
            smooth_x[nq+nλ+1:nq+nλ+nμ]    .= 0.0
            stepk! = stepk_maker(nq,nλ,nμ,qᵏ⁻¹,q̇ᵏ⁻¹,pᵏ⁻¹,tᵏ⁻¹,dynfuncs,invM,dt)

            for k = 1:iterations
                    stepk!(smooth_R,smooth_J,smooth_x,q̇ᵏ)
                    res = norm(smooth_R)
                    if res < ftol
                        isconverged = true
                        k_break = k-1
                        break
                    end
                    smooth_Δx = -smooth_J\smooth_R
                    smooth_x .+= smooth_Δx
                    qᵏ .= smooth_x[      1:nq]
                    λᵏ .= smooth_x[   nq+1:nq+nλ]
                    μᵏ .= smooth_x[nq+nλ+1:nq+nλ+nμ]
                    pᵏ .= Momentum_k(qᵏ⁻¹,pᵏ⁻¹,qᵏ,λᵏ,μᵏ,M,A,B,dt)
                    q̇ᵏ .= invM*pᵏ
            end
            q̃ᵏ .= smooth_x[      1:nq]
        else
            initial_x[            1:nq]             .= qᵏ⁻¹
            initial_x[         nq+1:nq+nq]          .= qᵏ⁻¹
            initial_x[      nq+nq+1:nq+nq+nλ]       .= 0.0
            initial_x[   nq+nq+nλ+1:nq+nq+nλ+nμ]    .= 0.0
            initial_x[nq+nq+nλ+nμ+1:nq+nq+nλ+nμ+nu] .= 0.0
            ns_stepk! = ns_stepk_maker(nq,nλ,nμ,nu,qᵏ⁻¹,q̇ᵏ⁻¹,pᵏ⁻¹,tᵏ⁻¹,Aset,dynfuncs,invM,dt)

            for k = 1:iterations

                nB, Bindex, B̄index = ns_stepk!(initial_R,initial_J,initial_x,q̇ᵏ)
                res = norm(initial_R)

                if res < ftol
                    isconverged = true
                    k_break = k-1
                    break
                else
                    Δx = -initial_J\initial_R
                    initial_x[1:nq+nq+nλ+nμ] .+= Δx[1:nq+nq+nλ+nμ]
                    initial_Λ = @view initial_x[nq+nq+nλ+nμ+1:nq+nq+nλ+nμ+nu]
                    ΔΛ = @view Δx[nq+nq+nλ+nμ+1:nq+nq+nλ+nμ+nu]
                    initial_Λ[Bindex] .+= ΔΛ[1:nB]
                    initial_Λ[B̄index] .+= ΔΛ[nB+1:nu]

                    qᵏ .= initial_x[         nq+1:nq+nq]
                    λᵏ .= initial_x[      nq+nq+1:nq+nq+nλ]
                    μᵏ .= initial_x[   nq+nq+nλ+1:nq+nq+nλ+nμ]
                    q̇ᵏ .= invM*Momentum_k(qᵏ⁻¹,pᵏ⁻¹,qᵏ,λᵏ,μᵏ,M,A,B,dt)
                end

                if count(Aset) != 0
                    @show k,res
                end
            end
            q̃ᵏ .= initial_x[            1:nq]
            qᵏ .= initial_x[         nq+1:nq+nq]
            λᵏ .= initial_x[      nq+nq+1:nq+nq+nλ]
            μᵏ .= initial_x[   nq+nq+nλ+1:nq+nq+nλ+nμ]
            Λᵏ .= initial_x[nq+nq+nλ+nμ+1:nq+nq+nλ+nμ+nu]
            pᵏ .= Momentum_k(qᵏ⁻¹,pᵏ⁻¹,qᵏ,λᵏ,μᵏ,M,A,B,dt)
            q̇ᵏ .= invM*pᵏ
        end

        if !isconverged
            if exception
                error("NLsolve max iterations $iterations, at timestep=$timestep, Res=$(res)")
            else
                # intor.convergence = false
                break
            end
        else
            # @info "timestep=$timestep, iterations=$k_break"
        end
        #---------Step k finisher-----------
        step += 1
        #---------Step k finisher-----------
        if verbose
            @printf("Progress: %5.1f%%, step: %s, time: %s, iterations: %s \n", (
            ts[timestep+1]/totaltime*100), timestep+1, ts[timestep+1], R_stepk_result.iterations)
        end
        next!(prog)
    end
    ts,qs,q̇s,ps,λs,friction_coefficients
end

function snhsolve(prob,nq,nλ,nμ,nu,q0,q̇0;tspan,dt=0.01,ftol=1e-14,verbose=false,iterations=50,
                progress=true,exception=true)
    @unpack bot,dynfuncs = prob
    M,Φ,A,Ψ,B,F!,jacobians,contact_funcs = dynfuncs
    Jac_F!,Ψq,∂Aᵀλ∂q,∂Bᵀμ∂q = jacobians
    E,𝐠,𝐠𝐪,∂𝐠𝐪ᵀΛ∂q,∂𝐠𝐪q̇∂q = contact_funcs
    totaltime = tspan[end] - tspan[begin]
    totalstep = ceil(Int,totaltime/dt)
    ts = [tspan[begin]+(i-1)*dt for i in 1:totalstep+1]
    qs = [copy(q0) for i in 1:totalstep+1]
    q̇s = [copy(q̇0) for i in 1:totalstep+1]
    ps = [M*copy(q̇0) for i in 1:totalstep+1]
    λs = [zeros(eltype(q0),nλ) for i in 1:totalstep+1]
    friction_coefficients = [zeros(eltype(q0),nμ) for i in 1:totalstep+1]
    Λs = [zeros(eltype(q0),nu) for i in 1:totalstep+1]
    invM = inv(M)
    # F⁺ = zero(q0)
    step = 0
    nx = nq + nλ + nμ + nu
    initial_x = zeros(eltype(q0),nx)
    initial_R = zeros(eltype(q0),nx)
    initial_J = zeros(eltype(q0),nx,nx)
    mr = norm(M,Inf)
    scaling = 1

    iteration = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = 1:totalstep
        #---------Step k Control-----------
        # control!(intor,cache)
        #---------Step k Control-----------
        qᵏ⁻¹ = qs[timestep]
        q̇ᵏ⁻¹ = q̇s[timestep]
        pᵏ⁻¹ = ps[timestep]
        λᵏ⁻¹ = λs[timestep]
        μᵏ⁻¹ = friction_coefficients[timestep]
        Λᵏ⁻¹ = Λs[timestep]
        tᵏ⁻¹ = ts[timestep]
        qᵏ   = qs[timestep+1]
        q̇ᵏ   = q̇s[timestep+1]
        pᵏ   = ps[timestep+1]
        λᵏ   = λs[timestep+1]
        μᵏ   = friction_coefficients[timestep+1]
        Λᵏ   = Λs[timestep+1]
        q̇ᵏ .= q̇ᵏ⁻¹

        qˣ = qᵏ⁻¹ .+ dt./2 .*q̇ᵏ⁻¹
        Aset = 𝐠(qˣ) .< 0
        Āset = .!Aset
        isconverged = false
        res = typemax(eltype(qᵏ))
        k_break = 0
        initial_x[         1:nq]          .= qᵏ⁻¹
        initial_x[      nq+1:nq+nλ]       .= 0.0
        initial_x[   nq+nλ+1:nq+nλ+nμ]    .= 0.0
        initial_x[nq+nλ+nμ+1:nq+nλ+nμ+nu] .= 0.0
        ns_stepk! = sns_stepk_maker(nq,nλ,nμ,nu,qᵏ⁻¹,q̇ᵏ⁻¹,pᵏ⁻¹,tᵏ⁻¹,Aset,dynfuncs,invM,dt)

        for k = 1:iterations

            nB, Bindex, B̄index = ns_stepk!(initial_R,initial_J,initial_x,q̇ᵏ)
            res = norm(initial_R)

            if res < ftol
                isconverged = true
                k_break = k-1
                break
            else
                Δx = -initial_J\initial_R
                initial_x[1:nq+nλ+nμ] .+= Δx[1:nq+nλ+nμ]
                initial_Λ = @view initial_x[nq+nλ+nμ+1:nq+nλ+nμ+nu]
                ΔΛ = @view Δx[nq+nλ+nμ+1:nq+nλ+nμ+nu]
                initial_Λ[Bindex] .+= ΔΛ[   1:nB]
                initial_Λ[B̄index] .+= ΔΛ[nB+1:nu]

                qᵏ .= initial_x[      1:nq]
                λᵏ .= initial_x[   nq+1:nq+nλ]
                μᵏ .= initial_x[nq+nλ+1:nq+nλ+nμ]
                q̇ᵏ .= invM*Momentum_k(qᵏ⁻¹,pᵏ⁻¹,qᵏ,λᵏ,μᵏ,M,A,B,dt)
            end

            if count(Aset) != 0
                @show k,res
            end
        end
        qᵏ .= initial_x[         1:nq]
        λᵏ .= initial_x[      nq+1:nq+nλ]
        μᵏ .= initial_x[   nq+nλ+1:nq+nλ+nμ]
        Λᵏ .= initial_x[nq+nλ+nμ+1:nq+nλ+nμ+nu]
        pᵏ .= Momentum_k(qᵏ⁻¹,pᵏ⁻¹,qᵏ,λᵏ,μᵏ,M,A,B,dt)
        q̇ᵏ .= invM*pᵏ

        if !isconverged
            if exception
                error("NLsolve max iterations $iterations, at timestep=$timestep, Res=$(res)")
            else
                # intor.convergence = false
                break
            end
        else
            @info "timestep=$timestep, iterations=$k_break"
        end
        #---------Step k finisher-----------
        step += 1
        #---------Step k finisher-----------
        if verbose
            @printf("Progress: %5.1f%%, step: %s, time: %s, iterations: %s \n", (
            ts[timestep+1]/totaltime*100), timestep+1, ts[timestep+1], R_stepk_result.iterations)
        end
        next!(prog)
    end
    ts,qs,q̇s,ps,λs,friction_coefficients
end

function ipsolve(prob,nq,nλ,nμ,q0,q̇0;dt=0.01,ftol=1e-14,verbose=false,iterations=50,
                progress=true,exception=true)
    @unpack bot,tspan,dynfuncs,control!,restart = prob
    M,Φ,A,Ψ,B,F!,jacobians,contact_funcs = dynfuncs
    Jac_F!,Ψq,∂Aᵀλ∂q,∂Bᵀμ∂q = jacobians
    E,𝐠,𝐠𝐪,∂𝐠𝐪ᵀΛ∂q,∂𝐠𝐪q̇∂q = contact_funcs
    totaltime = tspan[end] - tspan[begin]
    totalstep = ceil(Int,totaltime/dt)
    ts = [tspan[begin]+(i-1)*dt for i in 1:totalstep+1]
    q̃s = [copy(q0) for i in 1:totalstep+1]
    qs = [copy(q0) for i in 1:totalstep+1]
    q̇s = [copy(q̇0) for i in 1:totalstep+1]
    ps = [M*copy(q̇0) for i in 1:totalstep+1]
    λs = [zeros(eltype(q0),nλ) for i in 1:totalstep+1]
    friction_coefficients = [zeros(eltype(q0),nμ) for i in 1:totalstep+1]
    invM = inv(M)
    # F⁺ = zero(q0)
    step = 0
    smooth_nx = nq + nλ + nμ
    smooth_x = zeros(eltype(q0),smooth_nx)
    smooth_R = zeros(eltype(q0),smooth_nx)
    smooth_J = zeros(eltype(q0),smooth_nx,smooth_nx)
    mr = norm(M,Inf)
    scaling = 1

    iteration = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = 1:totalstep
        #---------Step k Control-----------
        # control!(intor,cache)
        #---------Step k Control-----------
        q̃ᵏ⁻¹ = q̃s[timestep]
        qᵏ⁻¹ = qs[timestep]
        q̇ᵏ⁻¹ = q̇s[timestep]
        pᵏ⁻¹ = ps[timestep]
        λᵏ⁻¹ = λs[timestep]
        μᵏ⁻¹ = friction_coefficients[timestep]
        tᵏ⁻¹ = ts[timestep]
        q̃ᵏ   = q̃s[timestep+1]
        qᵏ   = qs[timestep+1]
        q̇ᵏ   = q̇s[timestep+1]
        pᵏ   = ps[timestep+1]
        λᵏ   = λs[timestep+1]
        μᵏ   = friction_coefficients[timestep+1]
        q̇ᵏ .= q̇ᵏ⁻¹

        qˣ = qᵏ⁻¹ .+ dt./2 .*q̇ᵏ⁻¹
        Aset = 𝐠(qˣ) .< 0
        Āset = .!Aset
        isconverged = false
        res = typemax(eltype(qᵏ))
        k_break = 0
        nu = count(Aset)
        if nu == 0
            smooth_x[      1:nq]          .= qᵏ⁻¹
            smooth_x[   nq+1:nq+nλ]       .= 0.0
            smooth_x[nq+nλ+1:nq+nλ+nμ]    .= 0.0
            stepk! = stepk_maker(nq,nλ,nμ,qᵏ⁻¹,q̇ᵏ⁻¹,pᵏ⁻¹,tᵏ⁻¹,dynfuncs,invM,dt)

            for k = 1:iterations
                stepk!(smooth_R,smooth_J,smooth_x,q̇ᵏ)
                res = norm(smooth_R)
                if res < ftol
                    isconverged = true
                    k_break = k-1
                    break
                end
                smooth_Δx = -smooth_J\smooth_R
                smooth_x .+= smooth_Δx
                qᵏ .= smooth_x[      1:nq]
                λᵏ .= smooth_x[   nq+1:nq+nλ]
                μᵏ .= smooth_x[nq+nλ+1:nq+nλ+nμ]
                pᵏ .= Momentum_k(qᵏ⁻¹,pᵏ⁻¹,qᵏ,λᵏ,μᵏ,M,A,B,dt)
                q̇ᵏ .= invM*pᵏ
            end
            q̃ᵏ .= smooth_x[      1:nq]
        else
            nx = nq + nq + nλ + nμ + 2nu
            initial_x = zeros(eltype(q0),nx)
            initial_R = zeros(eltype(q0),nx)
            initial_J = zeros(eltype(q0),nx,nx)

            initial_x[            1:nq]             .= qᵏ⁻¹
            initial_x[         nq+1:nq+nq]          .= qᵏ⁻¹
            initial_x[      nq+nq+1:nq+nq+nλ]       .= λᵏ⁻¹
            initial_x[   nq+nq+nλ+1:nq+nq+nλ+nμ]    .= μᵏ⁻¹
            initial_x[nq+nq+nλ+nμ+1:nq+nq+nλ+nμ+2nu] .= 1.0
            ip_ns_stepk! = ip_ns_stepk_maker(nq,nλ,nμ,nu,qᵏ⁻¹,q̇ᵏ⁻¹,pᵏ⁻¹,tᵏ⁻¹,Aset,dynfuncs,invM,dt)
            μ = 1.0
            for k = 1:iterations
                ip_ns_stepk!(initial_R,initial_J,initial_x,q̇ᵏ,μ)
                res = norm(initial_R)
                @show k,res
                @show initial_R
                if res < ftol
                    isconverged = true
                    k_break = k-1
                    break
                end
                qᵏ .= initial_x[         nq+1:nq+nq]
                λᵏ .= initial_x[      nq+nq+1:nq+nq+nλ]
                μᵏ .= initial_x[   nq+nq+nλ+1:nq+nq+nλ+nμ]
                pᵏ .= Momentum_k(qᵏ⁻¹,pᵏ⁻¹,qᵏ,λᵏ,μᵏ,M,A,B,dt)
                q̇ᵏ .= invM*pᵏ
            end
            q̃ᵏ .= initial_x[            1:nq]
        end

        if !isconverged
            if exception
                error("NLsolve max iterations $iterations, at timestep=$timestep, Res=$(res)")
            else
                # intor.convergence = false
                break
            end
        else
            @info "timestep=$timestep, iterations=$k_break"
        end
        #---------Step k finisher-----------
        step += 1
        #---------Step k finisher-----------
        if verbose
            @printf("Progress: %5.1f%%, step: %s, time: %s, iterations: %s \n", (
            ts[timestep+1]/totaltime*100), timestep+1, ts[timestep+1], R_stepk_result.iterations)
        end
        next!(prog)
    end
    ts,qs,q̇s,ps,λs,friction_coefficients
end
