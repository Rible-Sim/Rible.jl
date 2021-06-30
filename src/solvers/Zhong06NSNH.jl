
function Momentum_k(qᵏ⁻¹,pᵏ⁻¹,qᵏ,λᵏ,μᵏ,M,A,B,h)
    pᵏ = -pᵏ⁻¹ .+ 2/h.*M*(qᵏ.-qᵏ⁻¹) .+
        1/(2h).*(transpose(A(qᵏ))-transpose(A(qᵏ⁻¹)))*λᵏ .+
        1/(2h).*(transpose(B(qᵏ))-transpose(B(qᵏ⁻¹)))*μᵏ
end

function stepk_maker(nq,nλ,nμ,nu,qᵏ⁻¹,q̇ᵏ⁻¹,pᵏ⁻¹,tᵏ⁻¹,dyfuncs,invM,h)
    M,Φ,A,Ψ,B,F!,jacobians,contact_funcs = dyfuncs
    Jac_F!,Ψq,∂Aᵀλ∂q,∂Bᵀμ∂q = jacobians
    # E,𝐠,𝐠𝐪,∂𝐠𝐪ᵀΛ∂q,∂𝐠𝐪q̇∂q = contact_funcs
    n1 =  nq
    n2 = 2nq
    n3 = 2nq+nλ
    n4 = 2nq+nλ+nμ
    function stepk!(R,J,x)
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

        q̇ᵏ = invM*Momentum_k(qᵏ⁻¹,pᵏ⁻¹,qᵏ,λᵏ,μᵏ,M,A,B,h)

        Aᵏ⁻¹ = A(qᵏ⁻¹)
        Bᵏ⁻¹ = B(qᵏ⁻¹)
        Aᵏ = A(qᵏ)
        Bᵏ = B(qᵏ)
        ∂q̇ᵏ∂qᵏ = 2/h*I + 1/(2h).*invM*(∂Aᵀλ∂q(qᵏ,λᵏ) + ∂Bᵀμ∂q(qᵏ,μᵏ))
        ∂q̇ᵏ∂λᵏ = invM*transpose(Aᵏ-Aᵏ⁻¹)/(2h)
        ∂q̇ᵏ∂μᵏ = invM*transpose(Bᵏ-Bᵏ⁻¹)/(2h)

        R[   1:n1] .= -h.*pᵏ⁻¹ .+ M*(q̃ᵏ.-qᵏ⁻¹) .-
                        1/2 .*transpose(Aᵏ⁻¹)*λᵏ .-
                        1/2 .*transpose(Bᵏ⁻¹)*μᵏ .-
                        (h^2)/2 .*F⁺
        R[n1+1:n2] .= (2/h).*M*(qᵏ - q̃ᵏ)
        R[n2+1:n3] .= Φ(qᵏ)
        R[n3+1:n4] .= Ψ(qᵏ,q̇ᵏ)
        R[n4+1:n4+nu] .= Λᵏ

        J .= 0.0
        J[   1:n1,   1:n1] .=  M
        J[   1:n1,n1+1:n2] .= -h^2/2 .*(1/2 .*∂F∂q .+ 1/h.*∂F∂q̇)
        J[   1:n1,n2+1:n3] .= -1/2 .*transpose(Aᵏ⁻¹)
        J[   1:n1,n3+1:n4] .= -1/2 .*transpose(Bᵏ⁻¹)

        J[n1+1:n2,   1:n1] .= -(2/h).*M
        J[n1+1:n2,n1+1:n2] .=  (2/h).*M

        J[n2+1:n3,n1+1:n2] .=  Aᵏ

        J[n3+1:n4,n1+1:n2] .=  Ψq(qᵏ,q̇ᵏ)+Bᵏ*∂q̇ᵏ∂qᵏ
        J[n3+1:n4,n2+1:n3] .=  Bᵏ*∂q̇ᵏ∂λᵏ
        J[n3+1:n4,n3+1:n4] .=  Bᵏ*∂q̇ᵏ∂μᵏ

        J[n4+1:n4+nu,n4+1:n4+nu] .= Matrix(1I,nu,nu)
    end

    stepk!
end

function ns_stepk_maker(nq,nλ,nμ,nu,qᵏ⁻¹,q̇ᵏ⁻¹,pᵏ⁻¹,tᵏ⁻¹,Aset,dyfuncs,invM,h)
    M,Φ,A,Ψ,B,F!,jacobians,contact_funcs = dyfuncs
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



function nhsolve(prob,nq,nλ,nμ,nu,q0,q̇0;dt=0.01,ftol=1e-14,verbose=false,iterations=50,
                progress=true,exception=true)
    @unpack bot,tspan,dyfuncs,control!,restart = prob
    M,Φ,A,Ψ,B,F!,jacobians,contact_funcs = dyfuncs
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
    μs = [zeros(eltype(q0),nμ) for i in 1:totalstep+1]
    Λs = [zeros(eltype(q0),nu) for i in 1:totalstep+1]
    invM = inv(M)
    # F⁺ = zero(q0)
    step = 0
    nx = nq + nq + nλ + nμ + nu
    initial_x = zeros(nx)
    initial_R = similar(initial_x)
    initial_J = zeros(eltype(initial_x),nx,nx)
    mr = norm(M,Inf)
    scaling = 1

    # ∂F∂q = zeros(eltype(q0),nq,nq)
    # ∂F∂q̇ = zeros(eltype(q0),nq,nq)



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
        μᵏ⁻¹ = μs[timestep]
        Λᵏ⁻¹ = Λs[timestep]
        tᵏ⁻¹ = ts[timestep]
        q̃ᵏ   = q̃s[timestep+1]
        qᵏ   = qs[timestep+1]
        q̇ᵏ   = q̇s[timestep+1]
        pᵏ   = ps[timestep+1]
        λᵏ   = λs[timestep+1]
        μᵏ   = μs[timestep+1]
        Λᵏ   = Λs[timestep+1]
        q̇ᵏ .= q̇ᵏ⁻¹
        initial_x[            1:nq]             .= qᵏ⁻¹
        initial_x[         nq+1:nq+nq]          .= qᵏ⁻¹
        initial_x[      nq+nq+1:nq+nq+nλ]       .= 0.0
        initial_x[   nq+nq+nλ+1:nq+nq+nλ+nμ]    .= 0.0
        initial_x[nq+nq+nλ+nμ+1:nq+nq+nλ+nμ+nu] .= 0.0
        # initial_R = similar(initial_x)
        #R_stepk!(initial_R,initial_x)
        # @show initial_R
        #@code_warntype R_stepk!(initial_R,initial_x)
        qˣ = qᵏ⁻¹ .+ dt./2 .*q̇ᵏ⁻¹
        Aset = 𝐠(qˣ) .< 0
        Āset = .!Aset
        ns_stepk! = ns_stepk_maker(nq,nλ,nμ,nu,qᵏ⁻¹,q̇ᵏ⁻¹,pᵏ⁻¹,tᵏ⁻¹,Aset,dyfuncs,invM,dt)
        stepk! = stepk_maker(nq,nλ,nμ,nu,qᵏ⁻¹,q̇ᵏ⁻¹,pᵏ⁻¹,tᵏ⁻¹,dyfuncs,invM,dt)
        isconverged = false
        res = typemax(eltype(qᵏ))
        k_break = 0
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
                # @show initial_R
                #
                # q̃ᵏ .= initial_x[            1:nq]
                # qᵏ .= initial_x[         nq+1:nq+nq]
                # # # @show q̃ᵏ
                # # # @show qᵏ
                # λᵏ .= initial_x[      nq+nq+1:nq+nq+nλ]
                # μᵏ .= initial_x[   nq+nq+nλ+1:nq+nq+nλ+nμ]
                # q̇ᵏ .= invM*Momentum_k(qᵏ⁻¹,pᵏ⁻¹,qᵏ,λᵏ,μᵏ,M,A,B,dt)
                # Λᵏ .= initial_x[nq+nq+nλ+nμ+1:nq+nq+nλ+nμ+nu]
                # @show q̇ᵏ⁻¹
                # @show q̇ᵏ
                # @show qᵏ⁻¹
                # @show q̃ᵏ
                # @show qᵏ
                # @show Λᵏ
            end
        end

        # for k = 1:iterations
        #
        #     stepk!(initial_R,initial_J,initial_x)
        #     res = norm(initial_R)
        #     if res < ftol
        #         isconverged = true
        #         k_break = k
        #         break
        #     else
        #         Δx = -initial_J\initial_R
        #         initial_x .+= Δx
        #         # @show k, res
        #         # @show initial_J[nq+1:nq+nq+nλ+nμ,nq+1:nq+nq+nλ+nμ]
        #     end
        # end
        # if converged(R_stepk_result)
        #     isconverged = true
        #     k_break = R_stepk_result.iterations
        #     initial_x[nq+1:nq+nq+nλ+nμ] = R_stepk_result.zero
        # end

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
        # iteration += R_stepk_result.iterations
        # @show R_stepk_result.iterations
        xᵏ = copy(initial_x)
        q̃ᵏ .= xᵏ[            1:nq]
        qᵏ .= xᵏ[         nq+1:nq+nq]
        λᵏ .= xᵏ[      nq+nq+1:nq+nq+nλ]
        μᵏ .= xᵏ[   nq+nq+nλ+1:nq+nq+nλ+nμ]
        Λᵏ .= xᵏ[nq+nq+nλ+nμ+1:nq+nq+nλ+nμ+nu]
        pᵏ .= Momentum_k(qᵏ⁻¹,pᵏ⁻¹,qᵏ,λᵏ,μᵏ,M,A,B,dt)
        q̇ᵏ .= invM*pᵏ
        # @show q̇ᵏ⁻¹
        # @show q̇ᵏ
        # Aq = A(qᵏ)
        # invM = inv(M+transpose(Aq)*Aq)
        #---------Step k finisher-----------
        step += 1
        # @show step
        # state.tprev = state.t
        # state.qprev .= state.q
        # state.q̇prev .= state.q̇
        # state.t = ts[timestep+1]
        # state.q .= qᵏ
        # state.q̇ .= q̇ᵏ
        #---------Step k finisher-----------
        if verbose
            @printf("Progress: %5.1f%%, step: %s, time: %s, iterations: %s \n", (
            ts[timestep+1]/totaltime*100), timestep+1, ts[timestep+1], R_stepk_result.iterations)
        end
        next!(prog)
    end
    ts,qs,q̇s,ps,λs,μs
end
