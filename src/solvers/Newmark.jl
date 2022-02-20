struct NewmarkCache{T,qT,λT,solverType}
    totalstep::Int
    totaltime::T
    ts::Vector{T}
    qs::qT
    q̇s::qT
    q̈s::qT
    λs::λT
    solver::solverType
end

function generate_cache(solver::Newmark,intor;dt,kargs...)
    @unpack prob,state,nx,nq,nλ = intor
    @unpack bot,tspan,dyfuncs,restart = prob
    @unpack q,q̇ = state.current
    @unpack M,F! = dyfuncs
    totaltime = tspan[end] - tspan[begin]
    totalstep = ceil(Int,totaltime/dt)
    ts = [tspan[begin]+(i-1)*dt for i in 1:totalstep+1]
    qs = [copy(q) for i in 1:totalstep+1]
    q̇s = [copy(q̇) for i in 1:totalstep+1]
    q̈s = [zero(q̇) for i in 1:totalstep+1]
    F0 = copy(q)
    F!(F0,qs[1],q̇s[1],0.0)
    q̈s[1] .= M\F0
    λs = [zeros(eltype(q),nλ) for i in 1:totalstep+1]
    # invM = inv(M)
    NewmarkCache(totalstep,totaltime,ts,qs,q̇s,q̈s,λs,solver)
end

function solve!(intor::Integrator,cache::NewmarkCache;
                dt=0.01,ftol=1e-14,verbose=false,iterations=50,
                progress=true,exception=true)
    @unpack prob,state,nx,nq,nλ = intor
    @unpack bot,tspan,dyfuncs,control!,restart = prob
    @unpack current,lasttime = state
    @unpack totaltime,totalstep,ts,qs,q̇s,q̈s,λs,solver = cache
    @unpack γ,β = solver
    @unpack M,Φ,A,F!,Jac_F! = dyfuncs
    q0 = qs[begin]
    F⁺ = zero(q0)
    F⁻ = zero(q0)
    step = 0
    initial_x = zeros(nx)
    initial_F = similar(initial_x)
    mr = norm(M,Inf)
    scaling = mr
    dt2 = dt^2
    d1 = 1/dt
    d2 = 1/dt^2
    b1 = 1/β
    function R_stepk_maker(qᵏ⁻¹,q̇ᵏ⁻¹,q̈ᵏ⁻¹,q̇ᵏ,q̈ᵏ,M,Φ,A,F!,nq,nλ,tᵏ⁻¹,dt)

        # @inline @inbounds
        function inner_R_stepk!(R,x)
            qᵏ = @view x[   1:nq]
            λᵏ = @view x[nq+1:nq+nλ]
            @. q̈ᵏ = ((qᵏ-qᵏ⁻¹)*d2 - q̇ᵏ⁻¹*d1)*b1 + (1-0.5*b1)*q̈ᵏ⁻¹
            @. q̇ᵏ = q̇ᵏ⁻¹ + ((1-γ)*q̈ᵏ⁻¹ + γ*q̈ᵏ)*dt
            F!(F⁺,qᵏ,q̇ᵏ,tᵏ⁻¹+dt)
            R[   1:nq]    .= dt2.*M*q̈ᵏ .+ transpose(A(qᵏ))*λᵏ .- dt2.*F⁺
            R[nq+1:nq+nλ] .= scaling.*Φ(qᵏ)
        end
    end


    iteration = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = 1:totalstep
        #---------Step k Control-----------
        control!(intor,cache)
        #---------Step k Control-----------
        qᵏ⁻¹ = qs[timestep]
        q̇ᵏ⁻¹ = q̇s[timestep]
        q̈ᵏ⁻¹ = q̈s[timestep]
        λᵏ⁻¹ = λs[timestep]
        tᵏ⁻¹ = ts[timestep]
        qᵏ = qs[timestep+1]
        q̇ᵏ = q̇s[timestep+1]
        q̈ᵏ = q̈s[timestep+1]
        λᵏ = λs[timestep+1]
        initial_x[   1:nq]    .= qᵏ⁻¹
        initial_x[nq+1:nq+nλ] .= λᵏ⁻¹
        # initial_R = similar(initial_x)
        #R_stepk!(initial_R,initial_x)
        # @show qᵏ⁻¹
        #@code_warntype R_stepk!(initial_R,initial_x)
        R_stepk! = R_stepk_maker(qᵏ⁻¹,q̇ᵏ⁻¹,q̈ᵏ⁻¹,q̇ᵏ,q̈ᵏ,M,Φ,A,F!,nq,nλ,tᵏ⁻¹,dt)
        if typeof(Jac_F!) == Nothing
            dfk = OnceDifferentiable(R_stepk!,initial_x,initial_F)
        else
            error("Jacobian not implemented yet.")
        end
        R_stepk_result = nlsolve(dfk, initial_x; ftol, iterations, method=:newton)
        if converged(R_stepk_result) == false
            if exception
                error("Not Converged!")
            else
                @warn "Step $timestep did not converged!"
                # intor.convergence = false
                # break
            end
        end
        iteration += R_stepk_result.iterations
        xᵏ = R_stepk_result.zero
        qᵏ .= xᵏ[   1:nq]
        λᵏ .= xᵏ[nq+1:nq+nλ]

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
            # dtwidth = ceil(Int,-log10(dt))
            @printf("Progress: %5.1f%%, step: %s, time: %.7f, iterations: %s \n", (
            ts[timestep+1]/totaltime*100), timestep+1, ts[timestep+1], R_stepk_result.iterations)
        end
        next!(prog)
    end

    return intor,cache
end

struct FBNewmarkCache{T,qT,λT,sT,solverType}
    totalstep::Int
    totaltime::T
    ts::Vector{T}
    qs::qT
    q̇s::qT
    q̈s::qT
    λs::λT
    s̄s::sT
    solver::solverType
end

function generate_cache(solver::FBNewmark,intor;dt,kargs...)
    @unpack prob,state,nx,nq,nλ = intor
    @unpack bot,tspan,dyfuncs,restart = prob
    @unpack q,q̇,s̄ = state.current
    @unpack M,F! = dyfuncs
    totaltime = tspan[end] - tspan[begin]
    totalstep = ceil(Int,totaltime/dt)
    ts = [tspan[begin]+(i-1)*dt for i in 1:totalstep+1]
    # ts = collect(range(tspan[begin];length=totalstep+1,step=dt))
    qs = [copy(q) for i in 1:totalstep+1]
    q̇s = [copy(q̇) for i in 1:totalstep+1]
    q̈s = [zero(q̇) for i in 1:totalstep+1]
    s̄s = [copy(s̄) for i in 1:totalstep+1]
    F0 = copy(q)
    F!(F0,qs[1],q̇s[1],s̄s[1],0.0)
    q̈s[1] .= M\F0
    λs = [zeros(eltype(q),nλ) for i in 1:totalstep+1]
    FBNewmarkCache(totalstep,totaltime,ts,qs,q̇s,q̈s,λs,s̄s,solver)
end

function solve!(intor::Integrator,cache::FBNewmarkCache;
                dt=0.01,ftol=1e-14,verbose=false,iterations=50,
                progress=true,exception=true)
    @unpack prob,state,nx,nq,nλ = intor
    @unpack bot,tspan,dyfuncs,control!,restart = prob
    @unpack current,lasttime = state
    @unpack totaltime,totalstep,ts,qs,q̇s,q̈s,s̄s,λs,solver = cache
    @unpack γ,β = solver.newmark
    @unpack M,Φ,A,F!,Ψ,apply_acu!,Jac_F! = dyfuncs
    q0 = qs[begin]
    F⁺ = zero(q0)
    F⁻ = zero(q0)
    step = 0
    initial_x = zeros(nx)
    initial_F = similar(initial_x)
    mr = norm(M,Inf)
    scaling = mr
    dt2 = dt^2
    d1 = 1/dt
    d2 = 1/dt^2
    b1 = 1/β
    function R_stepk_maker(qᵏ⁻¹,q̇ᵏ⁻¹,q̈ᵏ⁻¹,q̇ᵏ,q̈ᵏ,M,Φ,A,F!,Ψ,nq,nλ,tᵏ⁻¹,dt)

        function inner_R_stepk!(R,x)
            qᵏ = @view x[   1:nq]
            λᵏ = @view x[nq+1:nq+nλ]
            sᵏ = @view x[nq+nλ+1:nx]
            @. q̈ᵏ = ((qᵏ-qᵏ⁻¹)*d2 - q̇ᵏ⁻¹*d1)*b1 + (1-0.5*b1)*q̈ᵏ⁻¹
            @. q̇ᵏ = q̇ᵏ⁻¹ + ((1-γ)*q̈ᵏ⁻¹ + γ*q̈ᵏ)*dt
            F!(F⁺,qᵏ,q̇ᵏ,sᵏ,tᵏ⁻¹+dt)
            R[   1:nq]    .= M*q̈ᵏ .+ transpose(A(qᵏ))*λᵏ .- F⁺
            R[nq+1:nq+nλ] .= Φ(qᵏ)
            R[nq+nλ+1:nx] .= Ψ(sᵏ)
        end
    end

    function J_stepk_maker(qᵏ⁻¹,q̇ᵏ⁻¹,q̈ᵏ⁻¹,q̇ᵏ,q̈ᵏ,M,A,Jac_F!,nq,nλ,tᵏ⁻¹,dt)
        @inline @inbounds function inner_J_stepk!(J,x)
            qᵏ = @view x[   1:nq]
            λᵏ = @view x[nq+1:nq+nλ]
            s̄ = @view x[nq+nλ+1:end]
            q = qᵏ
            q̇ = q̇ᵏ⁻¹ + ((1-γ)*q̈ᵏ⁻¹ + γ*q̈ᵏ)*dt

            t = tᵏ⁻¹+ dt
            ∂F∂q,∂F∂q̇,∂F∂s̄ = Jac_F!(q,q̇,s̄,t)
            ζ = build_ζ(bot.tg)
            n = length(ζ)
            ∂Aᵀλ∂q1 = ∂Aᵀλ∂q(bot.tg,λᵏ)
            ∂ζ∂q = build_∂ζ∂q(bot.tg)(q)
            ∂ζ∂s̄ = build_∂ζ∂s̄(bot.tg)
            ∂s̄∂s̄ = I(n)
            κ₁ = 1000; κ₂ = 1000
            coes = diagm([((ζ[i]/κ₁)^2 + (κ₂*s̄[i])^2)^(-1/2) for i in 1:n])
            coζ = coes*diagm([ζ[i]/κ₁^2 for i in 1:n]) - diagm([1/κ₁ for i in 1:n])
            cos̄ = coes*diagm([κ₂^2*s̄[i] for i in 1:n]) - diagm([κ₂ for i in 1:n])
            #@show size(b1*d2*M ), size(diagm(λᵏ)*A(qᵏ)), size(∂F∂q), size(γ*b1*d1*∂F∂q̇)
            J[   1:nq ,      1:nq    ] .=  b1*d2*M  + ∂Aᵀλ∂q1 - ∂F∂q - γ*b1*d1*∂F∂q̇
            J[   1:nq ,   nq+1:nq+nλ ] .= transpose(A(qᵏ))
            J[   1:nq ,nq+nλ+1:end   ] .= -∂F∂s̄
            J[nq+1:nq+nλ,   1:nq ] .=  A(qᵏ)
            J[nq+1:nq+nλ,nq+1:nq+nλ] .=  0.0
            J[nq+1:nq+nλ,nq+nλ+1:end] .=  0.0
            J[nq+nλ+1:end,1:nq] .= coζ*∂ζ∂q
            J[nq+nλ+1:end,nq+1:nq+nλ+1] .= 0.0
            J[nq+nλ+1:end,nq+nλ+1:end] .= coζ*∂ζ∂s̄ + cos̄*∂s̄∂s̄
        end
    end

    function show_data(qᵏ⁻¹,sᵏ,M,A,Jac_F!,nq,nλ,tᵏ⁻¹,dt,initial_x,R_stepk!,J_stepk!)
        nx = length(initial_x)
        output = zeros(nx,nx)
        FiniteDiff.finite_difference_jacobian!(output,R_stepk!,initial_x)
        J = zeros(Float64, nx, nx)
        J_stepk!(J, initial_x)
        #wait()
        cha = J - output
        @show max = findmax(cha)
        max = 0
        if max > .02
            write_data(initial_x,output,J,cha) 
            wait()
        end
    end
    
    function write_data(initial_x,output,J,cha) 
        XLSX.openxlsx("t1.xlsx", mode="rw") do xf
            s1 = xf["1"]
            s2 = xf["2"]
            s3 = xf["3"]
            for i in 1:length(initial_x)
                s1[i,1:length(initial_x)] = output[i,:]
                s2[i,1:length(initial_x)] = J[i,:]
                s3[i,1:length(initial_x)] = cha[i,:]
            end
        end 
        @show 1
    end

    iteration = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = 1:totalstep
        #---------Step k Control-----------
        control!(intor,cache)
        #---------Step k Control-----------
        qᵏ⁻¹ = qs[timestep]
        q̇ᵏ⁻¹ = q̇s[timestep]
        q̈ᵏ⁻¹ = q̈s[timestep]
        λᵏ⁻¹ = λs[timestep]
        tᵏ⁻¹ = ts[timestep]
        sᵏ⁻¹ = s̄s[timestep]
        qᵏ = qs[timestep+1]
        q̇ᵏ = q̇s[timestep+1]
        q̈ᵏ = q̈s[timestep+1]
        λᵏ = λs[timestep+1]
        sᵏ = s̄s[timestep+1]
        initial_x[   1:nq]    .= qᵏ⁻¹
        initial_x[nq+1:nq+nλ] .= λᵏ⁻¹
        initial_x[nq+nλ+1:nx] .= sᵏ⁻¹
        # initial_R = similar(initial_x)
        #R_stepk!(initial_R,initial_x)
        #@show qᵏ⁻¹
        #@code_warntype R_stepk!(initial_R,initial_x)
        R_stepk! = R_stepk_maker(qᵏ⁻¹,q̇ᵏ⁻¹,q̈ᵏ⁻¹,q̇ᵏ,q̈ᵏ,M,Φ,A,F!,Ψ,nq,nλ,tᵏ⁻¹,dt)
        if typeof(Jac_F!) == Nothing
            dfk = OnceDifferentiable(R_stepk!,initial_x,initial_F)
        else
            J_stepk! = J_stepk_maker(qᵏ⁻¹,q̇ᵏ⁻¹,q̈ᵏ⁻¹,q̇ᵏ,q̈ᵏ,M,A,Jac_F!,nq,nλ,tᵏ⁻¹,dt)
            dfk = OnceDifferentiable(R_stepk!,J_stepk!,initial_x,initial_F)
        end

        if apply_acu! != nothing
            apply_acu!(bot.tg, tᵏ⁻¹;dt)
        end
        #show_data(qᵏ⁻¹,sᵏ,M,A,Jac_F!,nq,nλ,tᵏ⁻¹,dt,initial_x,R_stepk!,J_stepk!)

        R_stepk_result = nlsolve(dfk, initial_x; ftol, iterations, method=:newton)
        if converged(R_stepk_result) == false
            # @show R_stepk_result
            if exception
                error("Not Converged!")
            else
                @warn "Step $timestep did not converged!"
                # intor.convergence = false
                # break
            end
        end
        iteration += R_stepk_result.iterations
        xᵏ = R_stepk_result.zero
        qᵏ .= xᵏ[   1:nq]
        λᵏ .= xᵏ[nq+1:nq+nλ]
        sᵏ .= xᵏ[nq+nλ+1:nx]

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
            # dtwidth = ceil(Int,-log10(dt))
            @printf("Progress: %5.1f%%, step: %s, time: %.7f, iterations: %s \n", (
            ts[timestep+1]/totaltime*100), timestep+1, ts[timestep+1], R_stepk_result.iterations)
        end
        next!(prog)
    end

    return intor,cache
end


struct SNNewmarkCache{T,qT,λT,sT,solverType}
    totalstep::Int
    totaltime::T
    ts::Vector{T}
    qs::qT
    q̇s::qT
    q̈s::qT
    λs::λT
    s̄s::sT
    solver::solverType
end

function generate_cache(solver::SNNewmark,intor;dt,kargs...)
    @unpack prob,state,nx,nq,nλ = intor
    @unpack bot,tspan,dyfuncs,restart = prob
    @unpack q,q̇,s̄ = state.current
    @unpack M,F! = dyfuncs
    totaltime = tspan[end] - tspan[begin]
    totalstep = ceil(Int,totaltime/dt)
    ts = [tspan[begin]+(i-1)*dt for i in 1:totalstep+1]
    # ts = collect(range(tspan[begin];length=totalstep+1,step=dt))
    qs = [copy(q) for i in 1:totalstep+1]
    q̇s = [copy(q̇) for i in 1:totalstep+1]
    q̈s = [zero(q̇) for i in 1:totalstep+1]
    s̄s = [copy(s̄) for i in 1:totalstep+1]
    F0 = copy(q)
    F!(F0,qs[1],q̇s[1],s̄s[1],0.0)
    q̈s[1] .= M\F0
    λs = [zeros(eltype(q),nλ) for i in 1:totalstep+1]
    SNNewmarkCache(totalstep,totaltime,ts,qs,q̇s,q̈s,λs,s̄s,solver)
end
# using XLSX
function solve!(intor::Integrator,cache::SNNewmarkCache;
                dt=0.01,ftol=1e-14,verbose=false,iterations=50,
                progress=true,exception=true)
    @unpack prob,state,nx,nq,nλ = intor
    @unpack bot,tspan,dyfuncs,control!,restart = prob
    @unpack current,lasttime = state
    @unpack totaltime,totalstep,ts,qs,q̇s,q̈s,s̄s,λs,solver = cache
    @unpack γ,β = solver.newmark
    @unpack M,Φ,A,F!,Ψ,apply_acu!,Jac_F! = dyfuncs
    q0 = qs[begin]
    F⁺ = zero(q0)
    F⁻ = zero(q0)
    step = 0
    initial_x = zeros(nx)
    initial_F = similar(initial_x)
    mr = norm(M,Inf)
    scaling = mr
    dt2 = dt^2
    d1 = 1/dt
    d2 = 1/dt^2
    b1 = 1/β
    function R_stepk_maker(qᵏ⁻¹,q̇ᵏ⁻¹,q̈ᵏ⁻¹,q̇ᵏ,q̈ᵏ,M,Φ,A,F!,Ψ,nq,nλ,tᵏ⁻¹,dt)

        # @inline @inbounds
        function inner_R_stepk!(R,x)
            qᵏ = @view x[   1:nq]
            λᵏ = @view x[nq+1:nq+nλ]
            sᵏ = @view x[nq+nλ+1:nx]
            @. q̈ᵏ = ((qᵏ-qᵏ⁻¹)*d2 - q̇ᵏ⁻¹*d1)*b1 + (1-0.5*b1)*q̈ᵏ⁻¹
            @. q̇ᵏ = q̇ᵏ⁻¹ + ((1-γ)*q̈ᵏ⁻¹ + γ*q̈ᵏ)*dt
            F!(F⁺,qᵏ,q̇ᵏ,sᵏ,tᵏ⁻¹+dt)
            R[   1:nq]    .= M*q̈ᵏ .+ transpose(A(qᵏ))*λᵏ .- F⁺
            R[nq+1:nq+nλ] .= Φ(qᵏ)
            #R[nq+1:nq+nλ] .= Φ(qᵏ)
            ζ = build_ζ(bot.tg)
            add = nq+nλ
            c = 1
            n = length(ζ)
            i = findall(x->x<=0, ζ-c.*sᵏ)
            a = setdiff(1:n,i)
            na = length(a); ni = length(i)
            R[(1:ni) .+ add] = ζ[i]
            R[(1:na) .+ (ni+add)] = sᵏ[a]
        end
    end

    function J_stepk_maker(qᵏ⁻¹,q̇ᵏ⁻¹,q̈ᵏ⁻¹,q̇ᵏ,q̈ᵏ,M,A,Jac_F!,nq,nλ,tᵏ⁻¹,dt)
        @inline @inbounds function inner_J_stepk!(J,x)
            h = dt
            qᵏ = @view x[   1:nq]
            λᵏ = @view x[nq+1:nq+nλ]
            s̄ = @view x[nq+nλ+1:end]
            q = qᵏ
            q̇ = q̇ᵏ⁻¹ + ((1-γ)*q̈ᵏ⁻¹ + γ*q̈ᵏ)*dt

            t = tᵏ⁻¹+h
            ∂F∂q,∂F∂q̇,∂F∂s̄ = Jac_F!(q,q̇,s̄,t)
            ζ = build_ζ(bot.tg)
            c = 1
            n = length(ζ)
            i = findall(x->x<=0, ζ-c.*s̄)
            a = setdiff(1:n,i)
            ∂ζ∂q = build_∂ζ∂q(bot.tg)(q)[i,:]
            ∂ζ∂s̄ = build_∂ζ∂s̄(bot.tg)[i,:]
            ∂s̄∂s̄ = I(length(s̄))[a,:]
            ∂Aᵀλ∂q1 = ∂Aᵀλ∂q(bot.tg,λᵏ)

            #@show size(b1*d2*M ), size(diagm(λᵏ)*A(qᵏ)), size(∂F∂q), size(γ*b1*d1*∂F∂q̇)
            J[   1:nq ,      1:nq    ] .=  b1*d2*M  + ∂Aᵀλ∂q1 - ∂F∂q - γ*b1*d1*∂F∂q̇
            J[   1:nq ,   nq+1:nq+nλ ] .= transpose(A(qᵏ))
            J[   1:nq ,nq+nλ+1:end   ] .= -∂F∂s̄
            J[nq+1:nq+nλ,   1:nq ] .=  A(qᵏ)
            J[nq+1:nq+nλ,nq+1:nq+nλ] .=  0.0
            J[nq+1:nq+nλ,nq+nλ+1:end] .=  0.0
            #@show size(∂ζ∂q),size(zeros(Float64,length(a),nq))
            J[nq+nλ+1:end,1:nq] .= vcat(∂ζ∂q,zeros(Float64,length(a),nq))
            J[nq+nλ+1:end,nq+1:nq+nλ+1] .= 0.0
            J[nq+nλ+1:end,nq+nλ+1:end] .= vcat(∂ζ∂s̄,∂s̄∂s̄)
        end
    end

    function show_data(qᵏ⁻¹,sᵏ,M,A,Jac_F!,nq,nλ,tᵏ⁻¹,dt,initial_x,R_stepk!,J_stepk!)
        nx = length(initial_x)
        output = zeros(nx,nx)
        FiniteDiff.finite_difference_jacobian!(output,R_stepk!,initial_x)
        J = zeros(Float64, nx, nx)
        J_stepk!(J, initial_x)
        #wait()
        cha = J - output
        @show max = findmax(cha)
        max = 10
        if max > .02
            write_data(initial_x,output,J,cha) 
            wait()
        end
    end
    
    function write_data(initial_x,output,J,cha) 
        XLSX.openxlsx("t1.xlsx", mode="rw") do xf
            s1 = xf["1"]
            s2 = xf["2"]
            s3 = xf["3"]
            for i in 1:length(initial_x)
                s1[i,1:length(initial_x)] = output[i,:]
                s2[i,1:length(initial_x)] = J[i,:]
                s3[i,1:length(initial_x)] = cha[i,:]
            end
        end 
        @show 1
    end

    iteration = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    for timestep = 1:totalstep
        #---------Step k Control-----------
        control!(intor,cache)
        #---------Step k Control-----------
        qᵏ⁻¹ = qs[timestep]
        q̇ᵏ⁻¹ = q̇s[timestep]
        q̈ᵏ⁻¹ = q̈s[timestep]
        λᵏ⁻¹ = λs[timestep]
        tᵏ⁻¹ = ts[timestep]
        sᵏ⁻¹ = s̄s[timestep]
        qᵏ = qs[timestep+1]
        q̇ᵏ = q̇s[timestep+1]
        q̈ᵏ = q̈s[timestep+1]
        λᵏ = λs[timestep+1]
        sᵏ = s̄s[timestep+1]
        initial_x[   1:nq]    .= qᵏ⁻¹
        initial_x[nq+1:nq+nλ] .= λᵏ⁻¹
        initial_x[nq+nλ+1:nx] .= sᵏ⁻¹
        # initial_R = similar(initial_x)
        #R_stepk!(initial_R,initial_x)
        #@show qᵏ⁻¹
        #@code_warntype R_stepk!(initial_R,initial_x)
        R_stepk! = R_stepk_maker(qᵏ⁻¹,q̇ᵏ⁻¹,q̈ᵏ⁻¹,q̇ᵏ,q̈ᵏ,M,Φ,A,F!,Ψ,nq,nλ,tᵏ⁻¹,dt)
        if typeof(Jac_F!) == Nothing
            dfk = OnceDifferentiable(R_stepk!,initial_x,initial_F)
        else
            J_stepk! = J_stepk_maker(qᵏ⁻¹,q̇ᵏ⁻¹,q̈ᵏ⁻¹,q̇ᵏ,q̈ᵏ,M,A,Jac_F!,nq,nλ,tᵏ⁻¹,dt)
            dfk = OnceDifferentiable(R_stepk!,J_stepk!,initial_x,initial_F)
        end

        if apply_acu! != nothing
            apply_acu!(bot.tg, tᵏ⁻¹;dt)
        end
        #show_data(qᵏ⁻¹,sᵏ,M,A,Jac_F!,nq,nλ,tᵏ⁻¹,dt,initial_x,R_stepk!,J_stepk!)

        R_stepk_result = nlsolve(dfk, initial_x; ftol, iterations, method=:newton)
        if converged(R_stepk_result) == false
            # @show R_stepk_result
            if exception
                error("Not Converged!")
            else
                @warn "Step $timestep did not converged!"
                # intor.convergence = false
                # break
            end
        end
        iteration += R_stepk_result.iterations
        xᵏ = R_stepk_result.zero
        qᵏ .= xᵏ[   1:nq]
        λᵏ .= xᵏ[nq+1:nq+nλ]
        sᵏ .= xᵏ[nq+nλ+1:nx]

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
            # dtwidth = ceil(Int,-log10(dt))
            @printf("Progress: %5.1f%%, step: %.0f, time: %.3f, iterations: %.0f \n", (
            ts[timestep+1]/totaltime*100), timestep+1, ts[timestep+1], R_stepk_result.iterations)
        end
        next!(prog)
    end

    return intor,cache
end

