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
        @show qᵏ⁻¹
        #@code_warntype R_stepk!(initial_R,initial_x)
        R_stepk! = R_stepk_maker(qᵏ⁻¹,q̇ᵏ⁻¹,q̈ᵏ⁻¹,q̇ᵏ,q̈ᵏ,M,Φ,A,F!,nq,nλ,tᵏ⁻¹,dt)
        if typeof(Jac_F!) == Nothing
            dfk = OnceDifferentiable(R_stepk!,initial_x,initial_F)
        else
            error("Jacobian not implemented yet.")
        end
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

struct SlidingNewmarkCache{T,qT,λT,sT,solverType}
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

function generate_cache(solver::SlidingNewmark,intor;dt,kargs...)
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
    SlidingNewmarkCache(totalstep,totaltime,ts,qs,q̇s,q̈s,λs,s̄s,solver)
end

function solve!(intor::Integrator,cache::SlidingNewmarkCache;
                dt=0.01,ftol=1e-14,verbose=false,iterations=50,
                progress=true,exception=true)
    @unpack prob,state,nx,nq,nλ = intor
    @unpack bot,tspan,dyfuncs,control!,restart = prob
    @unpack current,lasttime = state
    @unpack totaltime,totalstep,ts,qs,q̇s,q̈s,s̄s,λs,solver = cache
    @unpack γ,β = solver.newmark
    @unpack M,Φ,A,F!,Ψ,Jac_F! = dyfuncs
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
            R[   1:nq]    .= dt2.*M*q̈ᵏ .+ transpose(A(qᵏ))*λᵏ .- dt2.*F⁺
            #R[   1:nq]    .= transpose(A(qᵏ))*λᵏ .- dt2.*F⁺
            R[nq+1:nq+nλ] .= scaling.*Φ(qᵏ)
            R[nq+nλ+1:nx] .= Ψ(sᵏ)
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
            error("Jacobian not implemented yet.")
        end
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
