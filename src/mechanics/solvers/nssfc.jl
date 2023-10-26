
module NSSFC
import ..Rible as TR
using Parameters
using LinearAlgebra
using StaticArrays
using BlockDiagonals
using Printf
using ProgressMeter
using OffsetArrays
using FiniteDiff
using FiniteDifferences

function ip_ns_stepk_maker(nq,nλ,nμ,nu,qₛ₋₁,q̇ₛ₋₁,pₛ₋₁,tₛ₋₁,dyfuncs,invM,h)
    M,Φ,A,Ψ,B,F!,jacobians,contact_funcs = dyfuncs
    Jac_F!,Ψq,∂Aᵀλ∂q,∂Bᵀμ∂q = jacobians
    𝐠,get_indices,get_FCs,get_D = contact_funcs

    stepk! = stepk_maker(nq,nλ,nμ,qₛ₋₁,q̇ₛ₋₁,pₛ₋₁,tₛ₋₁,dyfuncs,invM,h)

    nΛ = 3nu
    n1 = nq
    n2 = n1 + nλ
    n3 = n2 + nμ
    n4 = n3 + nΛ
    n5 = n4 + nΛ
    nx = n5
    T = eltype(qₛ₋₁)
    e = [one(T),zero(T),zero(T)]
    J = Diagonal([one(T),-one(T),-one(T)])
    𝐞_split = [e for i = 1:nu]
    function ip_ns_stepk!(𝐫𝐞𝐬,𝐉,x,Dₛ,ηs,es,H,μ)
        # function inner_res(x)
            qₛ = @view x[   1:n1]
            λₛ = @view x[n1+1:n2]
            μₛ = @view x[n2+1:n3]
            Λₛ = @view x[n3+1:n4]
            yₛ = @view x[n4+1:n5]

            𝐉 .= 0.0
            vₛ₋₁ = q̇ₛ₋₁

            q = (qₛ.+qₛ₋₁)./2
            q̇ = (qₛ.-qₛ₋₁)./h
            t = tₛ₋₁+h/2

            F⁺ = zeros(eltype(qₛ),nq)
            F!(F⁺,q,q̇,t)
            ∂F∂q = zeros(eltype(qₛ),nq,nq)
            ∂F∂q̇ = zeros(eltype(qₛ),nq,nq)
            Jac_F!(∂F∂q,∂F∂q̇,q,q̇,t)

            Aₛ₋₁ = A(qₛ₋₁)
            Bₛ₋₁ = B(qₛ₋₁)
            Aₛ = A(qₛ)
            Bₛ = B(qₛ)

            pₛ = Momentum_k(qₛ₋₁,pₛ₋₁,qₛ,λₛ,μₛ,M,A,B,h)
            vₛ = invM*pₛ

            ∂vₛ∂qₛ = 2/h*I + 1/(2h).*invM*(∂Aᵀλ∂q(qₛ,λₛ) + ∂Bᵀμ∂q(qₛ,μₛ))
            ∂vₛ∂λₛ = invM*transpose(Aₛ-Aₛ₋₁)/(2h)
            ∂vₛ∂μₛ = invM*transpose(Bₛ-Bₛ₋₁)/(2h)

            ∂Dₛvₛ∂qₛ = zeros(eltype(x),nΛ,nq)
            ∂bₛ∂qₛ = zeros(eltype(x),nΛ,nq)
            ∂DₛᵀHΛₛ∂qₛ = zeros(eltype(x),nq,nq)

            𝐛 = zeros(eltype(x),nΛ)
            ∂𝐛∂𝐱 = @view 𝐉[n3+1:n4,1:nx]
            v́ₛ₋₁ = Dₛ*vₛ₋₁
            v́ₛ = Dₛ*vₛ
            for i = 1:nu
                is = 3(i-1)
                vⁱₛ₋₁ = v́ₛ₋₁[is+1:is+3]
                vⁱₛ   = v́ₛ[is+1:is+3]
                vₜⁱₛ = norm(vⁱₛ[2:3])
                # vₜⁱₛ = sqrt(vⁱₛ[2]^2+vⁱₛ[3]^2+eps(eltype(x)))
                vᵤⁱₛ = norm(vⁱₛ[2])
                vᵥⁱₛ = norm(vⁱₛ[3])
                vₙⁱₛ₋₁ = vⁱₛ₋₁[1]
                vₙⁱₛ = vⁱₛ[1]
                # @show vₜⁱₛ, vₙⁱₛ₋₁, vₙⁱₛ
                𝐛[is+1] = ηs[i]*vₜⁱₛ + es[i]*vₙⁱₛ₋₁
                D̃i = copy(Dₛ[is+1:is+3,:])
                # D̃i[1,:] .+= (vⁱₛ[2].*Dₛ[is+2,:].+vⁱₛ[3].*Dₛ[is+3,:])./vₜⁱₛ
                ∂𝐛∂𝐱[is+1:is+3,   1:n1] .= D̃i*∂vₛ∂qₛ
                ∂𝐛∂𝐱[is+1:is+3,n1+1:n2] .= D̃i*∂vₛ∂λₛ
                ∂𝐛∂𝐱[is+1:is+3,n2+1:n3] .= D̃i*∂vₛ∂μₛ
            end
            # @show "before",yₛ,Λₛ
            Λ_split = RB.split_by_lengths(Λₛ,3)
            y_split = RB.split_by_lengths(yₛ,3)
            Λ_cone = [transpose(Λi)*J*Λi for Λi in Λ_split]
            y_cone = [transpose(yi)*J*yi for yi in y_split]
            # @show Λ_cone
            # @show y_cone
            # W_blocks = NTScale.(Λ_split,y_split)
            # z_split = W_blocks.*Λ_split
            # z = reduce(vcat,z_split)
            # W = BlockDiagonal(W_blocks)
            # WᵀW = transpose(W)*W
            # @show z_split
            # @show z_split⊙z_split
            𝐫𝐞𝐬 = zeros(eltype(x),nx)
            𝐫𝐞𝐬[   1:n1] .= -h.*pₛ₋₁ .+ M*(qₛ.-qₛ₋₁) .-
                            (2).*transpose(Dₛ)*H*Λₛ
                            1/2 .*transpose(Aₛ₋₁)*λₛ .-
                            1/2 .*transpose(Bₛ₋₁)*μₛ .-
                            (h^2)/2 .*F⁺
            𝐫𝐞𝐬[n1+1:n2] .= Φ(qₛ)
            𝐫𝐞𝐬[n2+1:n3] .= Ψ(qₛ,vₛ)
            𝐫𝐞𝐬[n3+1:n4] .= Dₛ*vₛ + 𝐛 - yₛ
            𝐫𝐞𝐬[n4+1:n5] .= reduce(vcat,Λ_split⊙y_split)

            # res = norm(𝐫𝐞𝐬)
            # @show 𝐫𝐞𝐬
            𝐉[   1:n1,   1:n1] .=  M .- (2).*∂DₛᵀHΛₛ∂qₛ .-(h^2/2).*(1/2 .*∂F∂q .+ 1/h.*∂F∂q̇)
            𝐉[   1:n1,n1+1:n2] .= -1/2 .*transpose(Aₛ₋₁)
            𝐉[   1:n1,n2+1:n3] .= -1/2 .*transpose(Bₛ₋₁)
            𝐉[   1:n1,n3+1:n4] .= -(2).*transpose(Dₛ)*H

            𝐉[n1+1:n2,   1:n1] .=  Aₛ

            𝐉[n2+1:n3,   1:n1] .=  Ψq(qₛ,vₛ) .+ Bₛ*∂vₛ∂qₛ
            𝐉[n2+1:n3,n1+1:n2] .=  Bₛ*∂vₛ∂λₛ
            𝐉[n2+1:n3,n2+1:n3] .=  Bₛ*∂vₛ∂μₛ

            𝐉[n3+1:n4,n4+1:n5] .= -Matrix(1I,nΛ,nΛ)

            𝐉[n4+1:n5,n3+1:n4] .=  BlockDiagonal(mat.(y_split))
            𝐉[n4+1:n5,n4+1:n5] .=  BlockDiagonal(mat.(Λ_split))
            # @show 𝐉
            # 𝐫𝐞𝐬
        # end

        # res = inner_res(x)
        # J_finite = FiniteDiff.finite_difference_jacobian(inner_res,x,Val{:central})
        # # J_finite = FiniteDifferences.jacobian(central_fdm(4, 1), inner_res, x)[1]
        # display(𝐉[n4+1:n5,n1+1:n2])
        # display(J_finite[n4+1:n5,n1+1:n2])
        # display(𝐉[n4+1:n5,n1+1:n2] - J_finite[n4+1:n5,n1+1:n2])

        # @show Λₛ,yₛ
        # η = 1.0
        Δxp = 𝐉\(-𝐫𝐞𝐬)
        ΔΛp = @view Δxp[n3+1:n4]
        Δyp = @view Δxp[n4+1:n5]
        ΔΛp_split = RB.split_by_lengths(ΔΛp,3)
        Δyp_split = RB.split_by_lengths(Δyp,3)
        # @show ΔΛp, Δyp
        # @show z_split,W_blocks,Δyp_split,ΔΛp_split,J
        αp_Λ = find_cone_step_length(Λ_split,ΔΛp_split,J)
        αp_y = find_cone_step_length(y_split,Δyp_split,J)
        αpmax = min(αp_Λ,αp_y)
        # αpmax = find_cone_step_length(z_split,W_blocks,Δyp_split,ΔΛp_split,J)
        αp = min(one(αpmax),0.99αpmax)
        Λp_split = Λ_split .+ αp.*ΔΛp_split
        yp_split = y_split .+ αp.*Δyp_split
        Λp_cone = [transpose(Λi)*J*Λi for Λi in Λp_split]
        yp_cone = [transpose(yi)*J*yi for yi in yp_split]
        # @show Λp_cone
        # @show yp_cone
        Λp = Λₛ .+ αp.*ΔΛp
        yp = yₛ .+ αp.*Δyp
        μp = transpose(yp)*Λp/nΛ
        σ = (μp/μ)^3
        τ = σ*μp
        # @show "Prediction",αpmax,αp,τ,σ,μ,μp
        # @show Δxp
        # @show αp.*ΔΛp,αp.*Δyp
        # @show Λp,yp
        𝐫𝐞𝐬_c_split = -σ.*μp.*𝐞_split.+((Δyp_split)⊙(ΔΛp_split))
        𝐫𝐞𝐬_c = reduce(vcat,𝐫𝐞𝐬_c_split)
        𝐫𝐞𝐬[n4+1:n5] .+= 𝐫𝐞𝐬_c
        Δxc = 𝐉\(-𝐫𝐞𝐬)
        # η = exp(-0.1μ) + 0.9
        ΔΛc = @view Δxc[n3+1:n4]
        Δyc = @view Δxc[n4+1:n5]
        ΔΛc_split = RB.split_by_lengths(ΔΛc,3)
        Δyc_split = RB.split_by_lengths(Δyc,3)
        # αmax = find_cone_step_length(z_split,W_blocks,Δyc_split,ΔΛc_split,J)
        α_Λ = find_cone_step_length(Λ_split,ΔΛc_split,J)
        # @show Λ_split,ΔΛc_split
        α_y = find_cone_step_length(y_split,Δyc_split,J)
        αmax = min(α_Λ,α_y)
        α = min(1,0.99αmax)
        qₛ .+= α.*Δxc[   1:n1]
        λₛ .+= α.*Δxc[n1+1:n2]
        μₛ .+= α.*Δxc[n2+1:n3]
        Λ_split .+= α.*ΔΛc_split
        y_split .+= α.*Δyc_split
        # @show ΔΛc_split, Δyc_split
        Λ_cone = [transpose(Λi)*J*Λi for Λi in Λ_split]
        y_cone = [transpose(yi)*J*yi for yi in y_split]
        Λₛ .+= α.*ΔΛc
        yₛ .+= α.*Δyc
        # @show "after",yₛ,Λₛ
        μ = transpose(yₛ)*Λₛ/nΛ
        # @show Λ_cone
        # @show y_cone
        # @show "Correction",αmax,α,τ,σ,μ,μp
        @show μ
        @show Λₛ
        # @show Δxc
        # @show α.*ΔΛc,α.*Δyc
        # @show Λₛ,yₛ
        μ,Δxc
    end
    ip_ns_stepk!
end

function ipsolve(nq,nλ,nμ,q0,q̇0,dyfuncs,tspan;dt=0.01,ftol=1e-14,xtol=ftol,verbose=false,imax=50,
                progress=true,exception=true)
    # @unpack bot,tspan,dyfuncs,control!,restart = prob
    M,Φ,A,Ψ,B,F!,jacobians,contact_funcs = dyfuncs
    Jac_F!,Ψq,∂Aᵀλ∂q,∂Bᵀμ∂q = jacobians
    𝐠,get_indices,get_FCs,get_D = contact_funcs
    totaltime = tspan[end] - tspan[begin]
    totalstep = ceil(Int,totaltime/dt)
    cs = OffsetArray([-1 for i in 1:totalstep+1],0:totalstep)
    ts = OffsetArray([tspan[begin]+(i-1)*dt for i in 1:totalstep+1],0:totalstep)
    q̃s = OffsetArray([copy(q0) for i in 1:totalstep+1],0:totalstep)
    qs = OffsetArray([copy(q0) for i in 1:totalstep+1],0:totalstep)
    q̇s = OffsetArray([copy(q̇0) for i in 1:totalstep+1],0:totalstep)
    ps = OffsetArray([M*copy(q̇0) for i in 1:totalstep+1],0:totalstep)
    λs = OffsetArray([zeros(eltype(q0),nλ) for i in 1:totalstep+1],0:totalstep)
    μs = OffsetArray([zeros(eltype(q0),nμ) for i in 1:totalstep+1],0:totalstep)
    # Λs = [zeros(eltype(q0),nu) for i in 1:totalstep+1]
    invM = inv(M)
    # F⁺ = zero(q0)
    step = 0
    smooth_nx = nq + nλ + nμ
    smooth_Δx = zeros(eltype(q0),smooth_nx)
    smooth_x = zero(smooth_Δx)
    smooth_R = zero(smooth_Δx)
    smooth_J = zeros(eltype(smooth_x),smooth_nx,smooth_nx)
    mr = norm(M,Inf)
    scaling = 1

    iteration = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = 1:totalstep
        #---------Step k Control-----------
        # control!(intor,cache)
        #---------Step k Control-----------
        q̃ₛ₋₁ = q̃s[timestep-1]
        qₛ₋₁ = qs[timestep-1]
        q̇ₛ₋₁ = q̇s[timestep-1]
        pₛ₋₁ = ps[timestep-1]
        λₛ₋₁ = λs[timestep-1]
        μₛ₋₁ = μs[timestep-1]
        tₛ₋₁ = ts[timestep-1]
        q̃ₛ   = q̃s[timestep]
        qₛ   = qs[timestep]
        q̇ₛ   = q̇s[timestep]
        pₛ   = ps[timestep]
        λₛ   = λs[timestep]
        μₛ   = μs[timestep]
        qˣ = qₛ₋₁ .+ dt./2 .*q̇ₛ₋₁
        qₛ .= qₛ₋₁ .+ dt.*q̇ₛ₋₁
        # q̇ₛ .= q̇ₛ₋₁
        nu,active_indices,g = get_indices(qˣ)
        gₙ = g[active_indices]
        cs[timestep] = nu
        # ns_stepk! = ns_stepk_maker(nq,nλ,nμ,qₛ₋₁,q̇ₛ₋₁,pₛ₋₁,tₛ₋₁,Aset,dyfuncs,invM,dt)
        # stepk! = stepk_maker(nq,nλ,nμ,nu,qₛ₋₁,q̇ₛ₋₁,pₛ₋₁,tₛ₋₁,dyfuncs,invM,dt)
        isconverged = false
        res = typemax(eltype(qₛ))
        iteration_break = 0
        if nu == 0
            smooth_x[      1:nq]          .= qₛ
            smooth_x[   nq+1:nq+nλ]       .= 0.0
            smooth_x[nq+nλ+1:nq+nλ+nμ]    .= 0.0
            stepk! = stepk_maker(nq,nλ,nμ,qₛ₋₁,q̇ₛ₋₁,pₛ₋₁,tₛ₋₁,dyfuncs,invM,dt)

            for iteration = 1:imax
                    stepk!(smooth_R,smooth_J,smooth_x)
                    res = norm(smooth_R)
                    if res < ftol
                        isconverged = true
                        iteration_break = iteration-1
                        break
                    end
                    smooth_Δx .= -smooth_J\smooth_R
                    smooth_x .+= smooth_Δx
            end
            qₛ .= smooth_x[      1:nq]
            λₛ .= smooth_x[   nq+1:nq+nλ]
            μₛ .= smooth_x[nq+nλ+1:nq+nλ+nμ]
            pₛ .= Momentum_k(qₛ₋₁,pₛ₋₁,qₛ,λₛ,μₛ,M,A,B,dt)
            q̇ₛ .= invM*pₛ
        else # u!=0
            nΛ = 3nu
            nonsmooth_nx = nq + nλ + nμ + nΛ + nΛ
            nonsmooth_x = zeros(eltype(q0),nonsmooth_nx)
            nonsmooth_R = zeros(eltype(q0),nonsmooth_nx)
            nonsmooth_J = zeros(eltype(q0),nonsmooth_nx,nonsmooth_nx)
            nonsmooth_x[            1:nq]                .= qₛ
            nonsmooth_x[   nq+nλ+nμ+1:nq+nλ+nμ+nΛ]    .= repeat([1,0,0],nu)
            nonsmooth_x[nq+nλ+nμ+nΛ+1:nq+nλ+nμ+nΛ+nΛ] .= repeat([1,0,0],nu)
            # y = copy(nonsmooth_x[nq+nq+nλ+nμ+1:nq+nq+nλ+nμ+nΛ])
            # q̇ₛ .= [0.5142598661572809, 0.0, 1.400342517987585]
            isconverged = false
            # @show timestep, nu
            ip_ns_stepk! = ip_ns_stepk_maker(nq,nλ,nμ,nu,qₛ₋₁,q̇ₛ₋₁,pₛ₋₁,tₛ₋₁,dyfuncs,invM,dt)
            μ = 1.0
            for iteration = 1:imax
                Dₛ,ηs,es,H = get_D(active_indices,qₛ)
                # ηs .= 1
                _,_,g = get_indices(qₛ)
                gₙ = g[active_indices]
                # @show iteration,Dₛ,ηs,es,gₙ
                μ,Δxc = ip_ns_stepk!(nonsmooth_R,nonsmooth_J,
                            nonsmooth_x,Dₛ,ηs,es,H,μ)
                res = norm(Δxc)
                @show timestep, iteration, res
                iteration_break = iteration
                if  res < ftol
                    isconverged = true
                    break
                end
            end
            qₛ .= nonsmooth_x[      1:nq]
            λₛ .= nonsmooth_x[   nq+1:nq+nλ]
            μₛ .= nonsmooth_x[nq+nλ+1:nq+nλ+nμ]
            pₛ .= Momentum_k(qₛ₋₁,pₛ₋₁,qₛ,λₛ,μₛ,M,A,B,dt)
            q̇ₛ .= invM*pₛ
            # @show gₙ*9.81
            # @show 𝚲ₛ./dt
        end

        if !isconverged
            @warn "NLsolve max iterations $iteration_break, at timestep=$timestep, Res=$(res)"
            if exception
                error("Not converged!")
            else
                # intor.convergence = false
                # break
            end
        else
            if nu == 0
                stepstring = "Smooth"
                # @info "$stepstring timestep=$timestep, iterations=$iteration_break"
            else
                stepstring = "Nonsmooth with $nu contact(s)"
                # @info "$stepstring timestep=$timestep, iterations=$iteration_break, and APGD_res=$APGD_res"
            end
        end

        #---------Step k finisher-----------
        step += 1
        #---------Step k finisher-----------
        if verbose
            @printf("Progress: %5.1f%%, step: %s, time: %s, iterations: %s \n", (
            ts[timestep]/totaltime*100), timestep, ts[timestep], R_stepk_result.iterations)
        end
        next!(prog)
    end
    ts,cs,qs,q̇s,ps,λs,μs
end

end
