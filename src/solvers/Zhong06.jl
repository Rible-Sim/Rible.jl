struct Zhong06Cache{T,qT,λT,MT}
    totalstep::Int
    totaltime::T
    ts::Vector{T}
    qs::qT
    q̇s::qT
    ps::qT
    λs::λT
    invM::MT
end

function generate_cache(::Zhong06,intor;dt,kargs...)
    @unpack prob,state,nx,nq,nλ = intor
    @unpack bot,tspan,dyfuncs,restart = prob
    @unpack current,lasttime = state
    @unpack q,q̇ = current
    @unpack M = dyfuncs
    totaltime = tspan[end] - tspan[begin]
    totalstep = ceil(Int,totaltime/dt)
    ts = [tspan[begin]+(i-1)*dt for i in 1:totalstep+1]
    qs = [copy(q) for i in 1:totalstep+1]
    q̇s = [copy(q̇) for i in 1:totalstep+1]
    ps = [M*copy(q̇) for i in 1:totalstep+1]
    λs = [zeros(eltype(q),nλ) for i in 1:totalstep+1]
    invM = inv(M)
    Zhong06Cache(totalstep,totaltime,ts,qs,q̇s,ps,λs,invM)
end

function solve!(intor::Integrator,cache::Zhong06Cache;
                dt=0.01,ftol=1e-14,verbose=false,iterations=50,
                progress=true,exception=true)
    @unpack prob,state,nx,nq,nλ = intor
    @unpack bot,tspan,dyfuncs,control!,restart = prob
    @unpack current,lasttime = state
    @unpack totaltime,totalstep,ts,qs,q̇s,ps,λs,invM = cache
    @unpack M,Φ,A,F!,Jac_F! = dyfuncs
    q0 = qs[begin]
    F⁺ = zero(q0)
    F⁻ = zero(q0)
    step = 0
    initial_x = zeros(nx)
    initial_F = similar(initial_x)
    mr = norm(M,Inf)
    scaling = mr

    @inline @inbounds function Momentum_k!(pᵏ,qᵏ⁻¹,pᵏ⁻¹,qᵏ,λᵏ,M,A,h)
        pᵏ .= -pᵏ⁻¹.+2/h.*M*(qᵏ.-qᵏ⁻¹) .+
            1/h.*scaling.*(transpose(A(qᵏ))-transpose(A(qᵏ⁻¹)))*λᵏ
    end
    function R_stepk_maker(qᵏ⁻¹,pᵏ⁻¹,M,Φ,A,F!,nq,nλ,tᵏ⁻¹,dt)
        # F!(F⁻,(qᵏ⁻¹.+qᵏ)./2,(qᵏ.-qᵏ⁻¹)./dt,tᵏ-dt/2)
        # kr = norm(F⁻,Inf)
        # scaling = mr + kr*dt^2
        @inline @inbounds function inner_R_stepk!(R,x)
            h = dt
            qᵏ = @view x[   1:nq]
            λᵏ = @view x[nq+1:nq+nλ]
            R[   1:nq]    .= M*(qᵏ.-qᵏ⁻¹) .-h.*pᵏ⁻¹.-scaling.*transpose(A(qᵏ⁻¹))*λᵏ
            F!(F⁺,(qᵏ.+qᵏ⁻¹)./2,(qᵏ.-qᵏ⁻¹)./h,tᵏ⁻¹+h/2)
            R[   1:nq]   .-= (h^2)/2 .*F⁺
            R[nq+1:nq+nλ] .= scaling.*Φ(qᵏ)
        end
    end

    ∂F∂q = zeros(eltype(q0),nq,nq)
    ∂F∂q̇ = zeros(eltype(q0),nq,nq)

    function J_stepk_maker(qᵏ⁻¹,M,A,Jac_F!,nq,nλ,tᵏ⁻¹,dt)
        @inline @inbounds function inner_J_stepk!(J,x)
            h = dt
            qᵏ = @view x[   1:nq]
            λᵏ = @view x[nq+1:end]
            q = (qᵏ.+qᵏ⁻¹)./2
            q̇ = (qᵏ.-qᵏ⁻¹)./h
            t = tᵏ⁻¹+h/2
            Jac_F!(∂F∂q,∂F∂q̇,q,q̇,t)
            J[   1:nq ,   1:nq]  .=  M.-h^2/2 .*(1/2 .*∂F∂q.+1/h.*∂F∂q̇)
            J[   1:nq ,nq+1:end] .= -scaling.*transpose(A(qᵏ⁻¹))
            J[nq+1:end,   1:nq ] .=  scaling.*A(qᵏ)
            J[nq+1:end,nq+1:end] .=  0.0
        end
    end

    iteration = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = 1:totalstep
        #---------Step k Control-----------
        control!(intor,cache)
        #---------Step k Control-----------
        qᵏ⁻¹ = qs[timestep]
        pᵏ⁻¹ = ps[timestep]
        λᵏ⁻¹ = λs[timestep]
        tᵏ⁻¹ = ts[timestep]
        qᵏ = qs[timestep+1]
        q̇ᵏ = q̇s[timestep+1]
        pᵏ = ps[timestep+1]
        λᵏ = λs[timestep+1]
        initial_x[   1:nq]    = qᵏ⁻¹
        initial_x[nq+1:nq+nλ] = λᵏ⁻¹
        # initial_R = similar(initial_x)
        #R_stepk!(initial_R,initial_x)
        # @show initial_R
        #@code_warntype R_stepk!(initial_R,initial_x)
        R_stepk! = R_stepk_maker(qᵏ⁻¹,pᵏ⁻¹,M,Φ,A,F!,nq,nλ,tᵏ⁻¹,dt)
        if typeof(Jac_F!) == Nothing
            dfk = OnceDifferentiable(R_stepk!,initial_x,initial_F)
        else
            J_stepk! = J_stepk_maker(qᵏ⁻¹,M,A,Jac_F!,nq,nλ,tᵏ⁻¹,dt)
            dfk = OnceDifferentiable(R_stepk!,J_stepk!,initial_x,initial_F)
        end
        R_stepk_result = nlsolve(dfk, initial_x; ftol, iterations, method=:newton)

        if converged(R_stepk_result) == false
            if exception
                error("Not Converged!")
            else
                # intor.convergence = false
                break
            end
        end
        iteration += R_stepk_result.iterations
        xᵏ = R_stepk_result.zero
        qᵏ .= xᵏ[   1:nq]
        λᵏ .= xᵏ[nq+1:nq+nλ]
        Momentum_k!(pᵏ,qᵏ⁻¹,pᵏ⁻¹,qᵏ,λᵏ,M,A,dt)

        # Aq = A(qᵏ)
        # invM = inv(M+transpose(Aq)*Aq)
        q̇ᵏ .= invM*pᵏ


        #---------Step k finisher-----------
        step += 1
        lasttime.t .= current.t
        lasttime.q .= current.q
        lasttime.q̇ .= current.q̇
        current.t .= ts[timestep+1]
        current.q .= qᵏ
        current.q̇ .= q̇ᵏ
        #---------Step k finisher-----------
        if verbose
            @printf("Progress: %5.1f%%, step: %s, time: %s, iterations: %s \n", (
            ts[timestep+1]/totaltime*100), timestep+1, ts[timestep+1], R_stepk_result.iterations)
        end
        next!(prog)
    end

    return intor,cache
end


struct FBZhong06Cache{T,qT,λT,sT,MT}
    totalstep::Int
    totaltime::T
    ts::Vector{T}
    qs::qT
    q̇s::qT
    ps::qT
    λs::λT
    s̄s::sT
    invM::MT
end

function generate_cache(::FBZhong06,intor;dt,kargs...)
    @unpack prob,state,nx,nq,nλ = intor
    @unpack bot,tspan,dyfuncs,restart = prob
    @unpack current,lasttime = state
    @unpack q,q̇,s̄ = current
    @unpack M = dyfuncs
    totaltime = tspan[end] - tspan[begin]
    totalstep = ceil(Int,totaltime/dt)
    ts = [tspan[begin]+(i-1)*dt for i in 1:totalstep+1]
    qs = [copy(q) for i in 1:totalstep+1]
    q̇s = [copy(q̇) for i in 1:totalstep+1]
    ps = [M*copy(q̇) for i in 1:totalstep+1]
    s̄s = [copy(s̄) for i in 1:totalstep+1]
    λs = [zeros(eltype(q),nλ) for i in 1:totalstep+1]
    invM = inv(M)
    FBZhong06Cache(totalstep,totaltime,ts,qs,q̇s,ps,λs,s̄s,invM)
end

function solve!(intor::Integrator,cache::FBZhong06Cache;
                dt=0.01,ftol=1e-14,verbose=false,iterations=50,
                progress=true,exception=true)
    @unpack prob,state,nx,nq,nλ = intor
    @unpack bot,tspan,dyfuncs,control!,restart = prob
    @unpack current,lasttime = state
    @unpack totaltime,totalstep,ts,qs,q̇s,ps,s̄s,λs,invM = cache
    @unpack M,Φ,A,F!,Ψ,apply_acu!,Jac_F! = dyfuncs
    q0 = qs[begin]
    F⁺ = zero(q0)
    F⁻ = zero(q0)
    step = 0
    initial_x = zeros(nx)
    initial_F = similar(initial_x)
    mr = norm(M,Inf)
    scaling = mr

    @inline @inbounds function Momentum_k!(pᵏ,qᵏ⁻¹,pᵏ⁻¹,qᵏ,λᵏ,M,A,h)
        pᵏ .= -pᵏ⁻¹.+2/h.*M*(qᵏ.-qᵏ⁻¹) .+
            1/h.*scaling.*(transpose(A(qᵏ))-transpose(A(qᵏ⁻¹)))*λᵏ
    end
    function R_stepk_maker(qᵏ⁻¹,pᵏ⁻¹,M,Φ,A,F!,nq,nλ,tᵏ⁻¹,dt)
        @inline @inbounds function inner_R_stepk!(R,x)
            h = dt
            qᵏ = @view x[   1:nq]
            λᵏ = @view x[nq+1:nq+nλ]
            sᵏ = @view x[nq+nλ+1:nx]
            R[   1:nq]    .= M*(qᵏ.-qᵏ⁻¹) .-h.*pᵏ⁻¹.-scaling.*transpose(A(qᵏ⁻¹))*λᵏ
            F!(F⁺,(qᵏ.+qᵏ⁻¹)./2,(qᵏ.-qᵏ⁻¹)./h,sᵏ,tᵏ⁻¹+h/2)
            #@show tᵏ⁻¹+h/2
            R[   1:nq]   .-= (h^2)/2 .*F⁺
            R[nq+1:nq+nλ] .= scaling.*Φ(qᵏ)
            R[nq+nλ+1:nx] .= Ψ(sᵏ)
        end
    end

    ∂F∂q = zeros(eltype(q0),nq,nq)
    ∂F∂q̇ = zeros(eltype(q0),nq,nq)

    function J_stepk_maker(qᵏ⁻¹,M,A,Jac_F!,nq,nλ,tᵏ⁻¹,dt)
        @inline @inbounds function inner_J_stepk!(J,x)
            h = dt
            qᵏ = @view x[   1:nq]
            λᵏ = @view x[nq+1:nq+nλ]
            s̄ = @view x[nq+nλ+1:end]
            q = (qᵏ.+qᵏ⁻¹)./2
            q̇ = (qᵏ.-qᵏ⁻¹)./h

            t = tᵏ⁻¹+h/2
            ∂F∂q,∂F∂q̇,∂F∂s̄ = Jac_F!(q,q̇,s̄,t)
            ζ = build_ζ(bot.tg)
            n = length(ζ)
            ∂ζ∂q = build_∂ζ∂q(bot.tg)(q)
            ∂ζ∂s̄ = build_∂ζ∂s̄(bot.tg)
            ∂s̄∂s̄ = I(n)
            κ₁ = 1000; κ₂ = 1000
            coes = diagm([((ζ[i]/κ₁)^2 + (κ₂*s̄[i])^2)^(-1/2) for i in 1:n])
            coζ = coes*diagm([ζ[i]/κ₁^2 for i in 1:n]) - diagm([1/κ₁ for i in 1:n])
            cos̄ = coes*diagm([κ₂^2*s̄[i] for i in 1:n]) - diagm([κ₂ for i in 1:n])
            J[   1:nq ,      1:nq    ] .=  M.-h^2/2 .*(1/2 .*∂F∂q.+1/h.*∂F∂q̇)
            J[   1:nq ,   nq+1:nq+nλ ] .= -scaling.*transpose(A(qᵏ⁻¹))
            J[   1:nq ,nq+nλ+1:end   ] .= -1/2 * dt^2 .* ∂F∂s̄
            J[nq+1:nq+nλ,   1:nq ] .=  scaling.*A(qᵏ)
            J[nq+1:nq+nλ,nq+1:nq+nλ] .=  0.0
            J[nq+1:nq+nλ,nq+nλ+1:end] .=  0.0
            #@show size(∂ζ∂q),size(zeros(Float64,length(a),nq))
            J[nq+nλ+1:end,1:nq] .= 1/2*coζ*∂ζ∂q
            J[nq+nλ+1:end,nq+1:nq+nλ+1] .= 0.0
            J[nq+nλ+1:end,nq+nλ+1:end] .= coζ *∂ζ∂s̄ + cos̄*∂s̄∂s̄
        end
    end


    iteration = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = 1:totalstep
        #---------Step k Control-----------
        control!(intor,cache)
        #---------Step k Control-----------
        qᵏ⁻¹ = qs[timestep]
        pᵏ⁻¹ = ps[timestep]
        λᵏ⁻¹ = λs[timestep]
        tᵏ⁻¹ = ts[timestep]
        sᵏ⁻¹ = s̄s[timestep]
        qᵏ = qs[timestep+1]
        q̇ᵏ = q̇s[timestep+1]
        pᵏ = ps[timestep+1]
        λᵏ = λs[timestep+1]
        sᵏ = s̄s[timestep+1]
        initial_x[   1:nq]    = qᵏ⁻¹
        initial_x[nq+1:nq+nλ] = λᵏ⁻¹
        initial_x[nq+nλ+1:nx] = sᵏ⁻¹

        if apply_acu! != nothing
            apply_acu!(bot.tg, tᵏ⁻¹;dt)
        end
        R_stepk! = R_stepk_maker(qᵏ⁻¹,pᵏ⁻¹,M,Φ,A,F!,nq,nλ,tᵏ⁻¹,dt)
        #dfk = OnceDifferentiable(R_stepk!,initial_x,initial_F)

        if typeof(Jac_F!) == Nothing
            dfk = OnceDifferentiable(R_stepk!,initial_x,initial_F)
        else
            J_stepk! = J_stepk_maker(qᵏ⁻¹,M,A,Jac_F!,nq,nλ,tᵏ⁻¹,dt)
            dfk = OnceDifferentiable(R_stepk!,J_stepk!,initial_x,initial_F)
        end

        #show_data(qᵏ⁻¹,sᵏ,M,A,Jac_F!,nq,nλ,tᵏ⁻¹,dt,initial_x,R_stepk!,J_stepk!)
        R_stepk_result = nlsolve(dfk, initial_x; ftol, iterations, method=:newton)
        #record_data(bot, R_stepk_result.iterations)

        if converged(R_stepk_result) == false
            if exception
                error("Not Converged!")
            else
                # intor.convergence = false
                @warn "Step $timestep did not converged!"
                #break
            end
        end
        iteration += R_stepk_result.iterations
        xᵏ = R_stepk_result.zero
        qᵏ .= xᵏ[   1:nq]
        λᵏ .= xᵏ[nq+1:nq+nλ]
        sᵏ .= xᵏ[nq+nλ+1:nx]
        Momentum_k!(pᵏ,qᵏ⁻¹,pᵏ⁻¹,qᵏ,λᵏ,M,A,dt)

        # Aq = A(qᵏ)
        # invM = inv(M+transpose(Aq)*Aq)
        q̇ᵏ .= invM*pᵏ


        #---------Step k finisher-----------
        step += 1
        lasttime.t .= current.t
        lasttime.q .= current.q
        lasttime.q̇ .= current.q̇
        lasttime.s̄ .= current.s̄
        current.t .= ts[timestep+1]
        current.q .= qᵏ
        current.q̇ .= q̇ᵏ
        current.s̄ .= sᵏ
        tension_lists = Vector{Vector{Float64}}()
        for clusterstring in bot.tg.clusterstrings
            tension_list = Vector{Float64}()
            for seg in clusterstring.segs
                push!(tension_list, seg.state.tension)
            end
            push!(tension_lists, tension_list)
        end
        push!(bot.traj.OtherData, tension_lists)

        #---------Step k finisher-----------
        if verbose
            @printf("Progress: %5.2f%%, step: %.0f, time: %.4f, iterations: %.0f \n", (
            ts[timestep+1]/totaltime*100), timestep+1, ts[timestep+1], R_stepk_result.iterations)
        end
        next!(prog)
    end

    return intor,cache
end

struct SNZhong06Cache{T,qT,λT,sT,MT}
    totalstep::Int
    totaltime::T
    ts::Vector{T}
    qs::qT
    q̇s::qT
    ps::qT
    λs::λT
    s̄s::sT
    invM::MT
end



function generate_cache(::SNZhong06,intor;dt,kargs...)
    @unpack prob,state,nx,nq,nλ = intor
    @unpack bot,tspan,dyfuncs,restart = prob
    @unpack current,lasttime = state
    @unpack q,q̇,s̄ = current
    @unpack M = dyfuncs
    totaltime = tspan[end] - tspan[begin]
    totalstep = ceil(Int,totaltime/dt)
    ts = [tspan[begin]+(i-1)*dt for i in 1:totalstep+1]
    qs = [copy(q) for i in 1:totalstep+1]
    q̇s = [copy(q̇) for i in 1:totalstep+1]
    ps = [M*copy(q̇) for i in 1:totalstep+1]
    s̄s = [copy(s̄) for i in 1:totalstep+1]
    λs = [zeros(eltype(q),nλ) for i in 1:totalstep+1]
    invM = inv(M)
    SNZhong06Cache(totalstep,totaltime,ts,qs,q̇s,ps,λs,s̄s,invM)
end

function solve!(intor::Integrator,cache::SNZhong06Cache;
                dt=0.01,ftol=1e-14,verbose=false,iterations=50,
                progress=true,exception=true)
    @unpack prob,state,nx,nq,nλ = intor
    @unpack bot,tspan,dyfuncs,control!,restart = prob
    @unpack current,lasttime = state
    @unpack totaltime,totalstep,ts,qs,q̇s,ps,s̄s,λs,invM = cache
    @unpack M,Φ,A,F!,Ψ,apply_acu!,Jac_F! = dyfuncs
    q0 = qs[begin]
    F⁺ = zero(q0)
    F⁻ = zero(q0)
    step = 0
    initial_x = zeros(nx)
    initial_F = similar(initial_x)
    mr = norm(M,Inf)
    scaling = mr

    @inline @inbounds function Momentum_k!(pᵏ,qᵏ⁻¹,pᵏ⁻¹,qᵏ,λᵏ,M,A,h)
        pᵏ .= -pᵏ⁻¹.+2/h.*M*(qᵏ.-qᵏ⁻¹) .+
            1/h.*scaling.*(transpose(A(qᵏ))-transpose(A(qᵏ⁻¹)))*λᵏ
    end
    function R_stepk_maker(qᵏ⁻¹,pᵏ⁻¹,M,Φ,A,F!,nq,nλ,tᵏ⁻¹,dt)
        # F!(F⁻,(qᵏ⁻¹.+qᵏ)./2,(qᵏ.-qᵏ⁻¹)./dt,tᵏ-dt/2)
        # kr = norm(F⁻,Inf)
        # scaling = mr + kr*dt^2
        @inline @inbounds function inner_R_stepk!(R,x)
            h = dt
            qᵏ = x[   1:nq]
            λᵏ = x[nq+1:nq+nλ]
            sᵏ = x[nq+nλ+1:nx]
            R[   1:nq]    .= M*(qᵏ.-qᵏ⁻¹) .-h.*pᵏ⁻¹.-scaling.*transpose(A(qᵏ⁻¹))*λᵏ
            F!(F⁺,(qᵏ.+qᵏ⁻¹)./2,(qᵏ.-qᵏ⁻¹)./h,sᵏ,tᵏ⁻¹+h/2)
            #@show tᵏ⁻¹+h/2
            R[   1:nq]   .-= (h^2)/2 .*F⁺
            R[nq+1:nq+nλ] .= scaling.*Φ(qᵏ)

            ζ = build_ζ(bot.tg)
            add = nq+nλ
            c = 1e8
            n = length(ζ)
            i = findall(x->x<=0, ζ-c.*sᵏ)
            a = setdiff(1:n,i)
            na = length(a); ni = length(i)
            R[(1:ni) .+ add] = ζ[i]
            R[(1:na) .+ (ni+add)] = sᵏ[a]
        end
    end

    ∂F∂q = zeros(eltype(q0),nq,nq)
    ∂F∂q̇ = zeros(eltype(q0),nq,nq)

    function J_stepk_maker(qᵏ⁻¹,M,A,Jac_F!,nq,nλ,tᵏ⁻¹,dt)
        @inline @inbounds function inner_J_stepk!(J,x)
            h = dt
            qᵏ = @view x[   1:nq]
            λᵏ = @view x[nq+1:nq+nλ]
            s̄ = @view x[nq+nλ+1:end]
            q = (qᵏ.+qᵏ⁻¹)./2
            q̇ = (qᵏ.-qᵏ⁻¹)./h

            t = tᵏ⁻¹+h/2
            ∂F∂q,∂F∂q̇,∂F∂s̄ = Jac_F!(q,q̇,s̄,t)
            ζ = build_ζ(bot.tg)
            c = 1e8
            n = length(ζ)
            i = findall(x->x<=0, ζ-c.*s̄)
            # @show ζ, s̄
            a = setdiff(1:n,i)
            ∂ζ∂q = build_∂ζ∂q(bot.tg)(q)[i,:]
            ∂ζ∂s̄ = build_∂ζ∂s̄(bot.tg)[i,:]
            ∂s̄∂s̄ = I(length(s̄))[a,:]
            J[   1:nq ,      1:nq    ] .=  M.-h^2/2 .*(1/2 .*∂F∂q.+1/h.*∂F∂q̇)
            J[   1:nq ,   nq+1:nq+nλ ] .= -scaling.*transpose(A(qᵏ⁻¹))
            J[   1:nq ,nq+nλ+1:end   ] .= -1/2 * dt^2 .* ∂F∂s̄
            J[nq+1:nq+nλ,   1:nq ] .=  scaling.*A(qᵏ)
            J[nq+1:nq+nλ,nq+1:nq+nλ] .=  0.0
            J[nq+1:nq+nλ,nq+nλ+1:end] .=  0.0
            #@show size(∂ζ∂q),size(zeros(Float64,length(a),nq))
            J[nq+nλ+1:end,1:nq] .= 1/2*vcat(∂ζ∂q,zeros(Float64,length(a),nq))
            J[nq+nλ+1:end,nq+1:nq+nλ+1] .= 0.0
            J[nq+nλ+1:end,nq+nλ+1:end] .= vcat(∂ζ∂s̄,∂s̄∂s̄)
            #@show cond(J)
        end
    end

    iteration = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = 1:totalstep
        #---------Step k Control-----------
        control!(intor,cache)
        #---------Step k Control-----------
        qᵏ⁻¹ = qs[timestep]
        pᵏ⁻¹ = ps[timestep]
        λᵏ⁻¹ = λs[timestep]
        tᵏ⁻¹ = ts[timestep]
        sᵏ⁻¹ = s̄s[timestep]
        q̇ᵏ⁻¹ = q̇s[timestep]
        qᵏ = qs[timestep+1]
        q̇ᵏ = q̇s[timestep+1]
        pᵏ = ps[timestep+1]
        λᵏ = λs[timestep+1]
        sᵏ = s̄s[timestep+1]
        initial_x[   1:nq]    = qᵏ⁻¹
        initial_x[nq+1:nq+nλ] = λᵏ⁻¹
        initial_x[nq+nλ+1:nx] = sᵏ⁻¹
        if apply_acu! != nothing
            apply_acu!(bot.tg, tᵏ⁻¹;dt)
        end
        # initial_R = similar(initial_x)
        #R_stepk!(initial_R,initial_x)
        # @show initial_R
        #@code_warntype R_stepk!(initial_R,initial_x)
        R_stepk! = R_stepk_maker(qᵏ⁻¹,pᵏ⁻¹,M,Φ,A,F!,nq,nλ,tᵏ⁻¹,dt)

        if typeof(Jac_F!) == Nothing
            dfk = OnceDifferentiable(R_stepk!,initial_x,initial_F)
        else
            J_stepk! = J_stepk_maker(qᵏ⁻¹,M,A,Jac_F!,nq,nλ,tᵏ⁻¹,dt)
            dfk = OnceDifferentiable(R_stepk!,J_stepk!,initial_x,initial_F)
        end

        R_stepk_result = nlsolve(dfk, initial_x; ftol, iterations, method=:newton)

        #xᵏ = newton(R_stepk!,J_stepk!,initial_x)
        #@show bot.tg.clusterstrings[1].segs[1].state.tension
        if converged(R_stepk_result) == false
            if exception
                error("Not Converged!")
            else
                # intor.convergence = false
                @warn "Step $timestep did not converged!"
                #break
            end
        end
        #iteration += R_stepk_result.iterations
        xᵏ = R_stepk_result.zero
        qᵏ .= xᵏ[   1:nq]
        λᵏ .= xᵏ[nq+1:nq+nλ]
        sᵏ .= xᵏ[nq+nλ+1:nx]
        Momentum_k!(pᵏ,qᵏ⁻¹,pᵏ⁻¹,qᵏ,λᵏ,M,A,dt)

        # Aq = A(qᵏ)
        # invM = inv(M+transpose(Aq)*Aq)
        q̇ᵏ .= invM*pᵏ


        #---------Step k finisher-----------
        step += 1
        lasttime.t .= current.t
        lasttime.q .= current.q
        lasttime.q̇ .= current.q̇
        lasttime.s̄ .= current.s̄
        current.t .= ts[timestep+1]
        current.q .= qᵏ
        current.q̇ .= q̇ᵏ
        current.s̄ .= sᵏ

        #---------Step k finisher-----------
        if verbose
            @printf("Progress: %5.1f%%, step: %.0f, time: %.3f, iterations: %.0f \n", (
            ts[timestep+1]/totaltime*100), timestep+1, ts[timestep+1], R_stepk_result.iterations)
        end
        #
        next!(prog)
    end

    return intor,cache
end

struct FBZhong06 end

struct FBZhong06Cache{MMT, AT, ΦT, ΨT}
    mass_matrices::MMT
    A::AT
    Φ::ΦT
    Ψ::ΨT
end

function generate_cache(::FBZhong06,intor;dt,kargs...)
    (;prob,state) = intor
    (;bot,dynfuncs) = prob
    (;q,q̇) = state.now
    # F!,_ = dynfuncs
    mm = build_MassMatrices(bot)
    # (;M) = mm
    A = make_A(bot)
    Φ = make_Φ(bot)
    Ψ = make_Ψ(bot.tg)
    FBZhong06Cache(mm,A,Φ,Ψ)
end

function solve!(intor::Integrator,cache::FBZhong06Cache;
                dt,ftol=1e-14,verbose=false,iterations=50,
                progress=true,exception=true)
    (;prob,state,control!,tspan,restart,totalstep) = intor
    (;bot,dynfuncs) = prob
    (;traj) = bot
    # @unpack t,q,q̇,tprev,qprev,q̇prev = state
    (;F!,Jac_F!,apply_acu!) = dynfuncs
    (;mass_matrices,A,Φ,Ψ) = cache
    (;Ḿ,M̌,M̄,invM̌) = mass_matrices
    q̌0 = traj.q̌[begin]
    λ0 = traj.λ[begin]
    q̇0 = traj.q̇[begin]
    s0 = traj.s[begin]
    T = eltype(q̌0)
    nq̌ = length(q̌0)
    nλ = length(λ0)
    ns = length(s0)
    ∂F∂q̌ = zeros(T,nq̌,nq̌)
    ∂F∂q̌̇ = zeros(T,nq̌,nq̌)
    p̌ᵏ⁻¹ = Ḿ*q̇0
    p̌ᵏ   = zero(p̌ᵏ⁻¹)
    step = 0
    initial_x = vcat(q̌0,λ0,s0)
    initial_Res = zero(initial_x)
    nx = length(initial_x)
    mr = norm(Ḿ,Inf)
    scaling = mr

    function make_Res_stepk(qᵏ,q̌ᵏ,λᵏ,sᵏ,qᵏ⁻¹,p̌ᵏ⁻¹,F̌,Aᵀ,tᵏ⁻¹)
        @inline @inbounds function inner_Res_stepk!(Res,x)
            h = dt
            q̌ᵏ .= x[   1:nq̌   ]
            λᵏ .= x[nq̌+1:nq̌+nλ]
            sᵏ .= x[nq̌+nλ+1:nx]
            F!(F̌,(qᵏ.+qᵏ⁻¹)./2,(qᵏ.-qᵏ⁻¹)./h,sᵏ,tᵏ⁻¹+h/2)
            Res[   1:nq̌   ] .= Ḿ*(qᵏ.-qᵏ⁻¹) .-
                               h.*p̌ᵏ⁻¹ .-
                               (h^2)/2 .*F̌ .+
                               scaling.*Aᵀ*λᵏ
            Res[nq̌+1:nq̌+nλ] .= scaling.*Φ(qᵏ)
            # Res[1:nq̌] .-= (h^2)/2 .* F̌
            Res[nq̌+nλ+1:nx] .= Ψ(sᵏ)
        end
    end



    function make_Jac_stepk(qᵏ,q̌ᵏ,sᵏ,qᵏ⁻¹,∂F∂q̌,∂F∂q̌̇,Aᵀ,tᵏ⁻¹)
        @inline @inbounds function inner_Jac_stepk!(Jac,x)
            h = dt
            q̌ᵏ .= x[1:nq̌]
            sᵏ .= x[nq̌+nλ+1:nx]
            ζ = build_ζ(bot.tg)
            ∂ζ∂q = build_∂ζ∂q(bot.tg, q̌ᵏ)
            ∂ζ∂s̄ = build_∂ζ∂s̄(bot.tg)
            n = length(ζ)
            ∂s̄∂s̄ = I(n)
            κ₁ = 10; κ₂ = 10
            coes = diagm([((ζ[i]/κ₁)^2 + (κ₂*sᵏ[i])^2)^(-1/2) for i in 1:n])
            coζ = coes*diagm([ζ[i]/κ₁^2 for i in 1:n]) - diagm([1/κ₁ for i in 1:n])
            cos̄ = coes*diagm([κ₂^2*sᵏ[i] for i in 1:n]) - diagm([κ₂ for i in 1:n])
            # Jac_F!(∂F∂q̌,∂F∂q̌̇,(qᵏ.+qᵏ⁻¹)./2,(qᵏ.-qᵏ⁻¹)./h,tᵏ⁻¹+h/2)
            ∂Q̌∂q̌̇ = build_∂Q̌∂q̌̇(bot.tg)
            ∂Q̌∂q̌ = build_∂Q̌∂q̌(bot.tg)
            ∂Q̌∂s̄ = build_∂Q̌∂s̄(bot.tg)
            Jac[   1:nq̌ ,      1:nq̌    ] .=  M̌.-h^2/2 .*(1/2 .*∂Q̌∂q̌.+1/h.*∂Q̌∂q̌̇)
            Jac[   1:nq̌ ,   nq̌+1:nq̌+nλ ] .= scaling.*Aᵀ
            # @show nq̌, length(Jac[1, nq̌+nλ+1:end])
            # @show size(∂Q̌∂s̄)
            Jac[   1:nq̌ ,   nq̌+nλ+1:end] .= -1/2 * dt^2 .* ∂Q̌∂s̄
            Jac[nq̌+1:nq̌+nλ,   1:nq̌ ] .=  scaling.*A(qᵏ)
            Jac[nq̌+1:nq̌+nλ,nq̌+1:nq̌+nλ] .=  0.0
            Jac[nq̌+1:nq̌+nλ,nq̌+nλ+1:end] .=  0.0
            Jac[nq̌+nλ+1:end,1:nq̌] .= 1/2 * coζ*∂ζ∂q
            Jac[nq̌+nλ+1:end,nq̌+1:nq̌+nλ+1] .= 0.0
            Jac[nq̌+nλ+1:end,nq̌+nλ+1:end] .= coζ *∂ζ∂s̄ + cos̄*∂s̄∂s̄
        end
    end

    function test(qᵏ,q̌ᵏ,sᵏ,qᵏ⁻¹,∂F∂q̌,∂F∂q̌̇,Aᵀ,tᵏ⁻¹)
        @inline @inbounds function inner_Jac_stepk!(x)
            h = dt
            q̌ᵏ .= x[1:nq̌]
            sᵏ .= x[nq̌+nλ+1:nx]
            ζ = build_ζ(bot.tg)
            ∂ζ∂q = Record_build_∂ζ∂q(bot.tg, q̌ᵏ,"t1.xlsx","t1")
            ∂ζ∂s̄ = build_∂ζ∂s̄(bot.tg)
            n = length(ζ)
            ∂s̄∂s̄ = I(n)
            κ₁ = 10; κ₂ = 10
            coes = diagm([((ζ[i]/κ₁)^2 + (κ₂*sᵏ[i])^2)^(-1/2) for i in 1:n])
            coζ = coes*diagm([ζ[i]/κ₁^2 for i in 1:n]) - diagm([1/κ₁ for i in 1:n])
            cos̄ = coes*diagm([κ₂^2*sᵏ[i] for i in 1:n]) - diagm([κ₂ for i in 1:n])
            # Jac_F!(∂F∂q̌,∂F∂q̌̇,(qᵏ.+qᵏ⁻¹)./2,(qᵏ.-qᵏ⁻¹)./h,tᵏ⁻¹+h/2)
            ∂Q̌∂q̌̇ = build_∂Q̌∂q̌̇(bot.tg)
            ∂Q̌∂q̌ = build_∂Q̌∂q̌(bot.tg)
            ∂Q̌∂s̄ = build_∂Q̌∂s̄(bot.tg)
            # @show nq̌, length(Jac[1, nq̌+nλ+1:end])
            # @show size(∂Q̌∂s̄)
            return coζ, ∂ζ∂q, q̌ᵏ
        end
    end

    @inline @inbounds function Momentum_k!(p̌ᵏ,p̌ᵏ⁻¹,qᵏ,qᵏ⁻¹,λᵏ,M,A,Aᵀ,h)
        p̌ᵏ .= -p̌ᵏ⁻¹.+2/h.*Ḿ*(qᵏ.-qᵏ⁻¹) .-
            1/h.*scaling.*(transpose(A(qᵏ))-Aᵀ)*λᵏ
    end

    total_iterations = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = 1:totalstep
        #---------Step k Control-----------
        # control!(intor,cache)
        #---------Step k Control-----------
        tᵏ⁻¹ = traj.t[timestep]
        tᵏ = traj.t[timestep+1]
        qᵏ⁻¹ = traj.q[timestep]
        sᵏ⁻¹ = traj.s[timestep]
        sᵏ = traj.s[timestep+1]
        qᵏ = traj.q[timestep+1]
        q̌ᵏ = traj.q̌[timestep+1]
        q̌̇ᵏ = traj.q̌̇[timestep+1]
        q̃̇ᵏ = traj.q̃̇[timestep+1]
        λᵏ = traj.λ[timestep+1]
        F̌ = traj.F̌[timestep+1]
        initial_x[   1:nq̌]    .= traj.q̌[timestep]
        initial_x[nq̌+1:nq̌+nλ] .= traj.λ[timestep]
        initial_x[nq̌+nλ+1:nx] .= traj.s[timestep]
        initial_x
        Aᵀ = transpose(A(qᵏ⁻¹))
        if apply_acu! != nothing
            apply_acu!(bot.tg, tᵏ⁻¹; dt)
        end
        Res_stepk! = make_Res_stepk(qᵏ,q̌ᵏ,λᵏ,sᵏ,qᵏ⁻¹,p̌ᵏ⁻¹,F̌,Aᵀ,tᵏ⁻¹)
        if Jac_F! isa Nothing
            dfk = OnceDifferentiable(Res_stepk!,initial_x,initial_Res)
        else

            Jac_stepk! = make_Jac_stepk(qᵏ,q̌ᵏ,sᵏ,qᵏ⁻¹,∂F∂q̌,∂F∂q̌̇,Aᵀ,tᵏ⁻¹)
            dfk = OnceDifferentiable(Res_stepk!,Jac_stepk!,initial_x,initial_Res)
        end

        # if timestep==totalstep
        #     nx = length(initial_x)
        #     output = zeros(nx, nx)
        #     FiniteDiff.finite_difference_jacobian!(output, Res_stepk!, initial_x)
        #     J = zeros(Float64, nx, nx)
        #     J[findall(x->abs(x)<1e-7, J)] .= 0
        #     Jac_stepk!(J, initial_x)
        #     diff = J - output

        #     test! = test(qᵏ,q̌ᵏ,sᵏ,qᵏ⁻¹,∂F∂q̌,∂F∂q̌̇,Aᵀ,tᵏ⁻¹)
        #     global coζ, dζdq, q̌ = test!(initial_x)
        #     h1,w1 = size(coζ); h2, w2 = size(dζdq)
        #     w3 = length(q̌)
        #     XLSX.openxlsx("t1.xlsx", mode="rw") do xf
        #         s1 = xf["前向差分"]
        #         s2 = xf["精确值"]
        #         s3 = xf["两者之差"]
        #         s4 = xf["Temp"]
        #         for i in 1:nx
        #             s1[i, 1:nx] = output[i, :]
        #             s2[i, 1:nx] = J[i, :]
        #             s3[i, 1:nx] = diff[i, :]
        #         end
        #         for i in 1:h1
        #             s4[i, 1:w1] = coζ[i, :]
        #         end
        #         for i in h1+1:h1+h2
        #             s4[i, 1:w2] = dζdq[i-h1, :]
        #         end
        #         for i in h1+h2+1
        #             s4[i, 1:w3] = q̌[1:w3]
        #         end
        #     end
        # end

        Res_stepk_result = nlsolve(dfk, initial_x; ftol, iterations, method=:newton)

        if converged(Res_stepk_result) == false
            if exception
                error("Not Converged! Step=$timestep")
            else
                # intor.convergence = false
                break
            end
        end
        total_iterations += Res_stepk_result.iterations
        xᵏ = Res_stepk_result.zero
        q̌ᵏ .= xᵏ[   1:nq̌   ]
        λᵏ .= xᵏ[nq̌+1:nq̌+nλ]
        sᵏ .= xᵏ[nq̌+nλ+1:nx]
        Momentum_k!(p̌ᵏ,p̌ᵏ⁻¹,qᵏ,qᵏ⁻¹,λᵏ,Ḿ,A,Aᵀ,dt)
        q̌̇ᵏ .= invM̌*(p̌ᵏ.-M̄*q̃̇ᵏ )
        #---------Step k finisher-----------
        step += 1
        p̌ᵏ,p̌ᵏ⁻¹ = p̌ᵏ⁻¹,p̌ᵏ
        state.prv = traj[timestep]
        state.now = traj[timestep+1]
        #---------Step k finisher-----------
        if verbose
            dg_step = ceil(Int,log10(totalstep))+1
            dg_dt = max(1,-floor(Int,log10(dt)))
            wd_t = ceil(Int,log10(traj.t[end]))+dg_dt+1+1
            progfmt = Printf.Format("Progress: %5.1f%%, step: %$(dg_step)u, time: %$(wd_t).$(dg_dt)f, iterations: %s \n")
            progstr = Printf.format(progfmt,
                floor(timestep/totalstep*100;digits=1), timestep, tᵏ, Res_stepk_result.iterations
            )
            print(progstr)
        end
        next!(prog)
    end

    return intor,cache
end
