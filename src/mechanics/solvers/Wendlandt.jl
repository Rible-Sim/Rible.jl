
struct WendlandtCache end
generate_cache(::Wendlandt,sim) = WendlandtCache()

function solve(sim::Simulator,cache::WendlandtCache;dt=0.01,ftol=1e-14,verbose=false,
                                 callback=DEFAULT_CALLBACK)
    (;prob,sol) = sim
    (;ts,qs,q̇s,λs) = sol
    (;funcs,tspan,q0,q̇0,λ0,nx,nq,nλ) = prob
    M,Φ,A,F!,Jacs = funcs
    F⁺ = zero(q0)
    F⁻ = zero(q0)
    tend = tspan[2]
    step = 0
    initial_x = zeros(nx)
    mr = norm(M,Inf)
    scaling = mr
    # First step
    function R_step1_maker(q0,q̇0,M,Φ,A,F!,nq,nλ,dt)
        # F!(F⁺,q0,q̇0,0.0)
        # kr = norm(F⁺,Inf)
        # scaling = mr + kr*dt^2
        function inner_R_step1!(R,x)
            h = dt
            myq1 = @view x[   1:nq]
            myλ0 = @view x[nq+1:nq+nλ]
            R[   1:nq]    .= M*((myq1.-q0) .- q̇0.*h) .- transpose(A(q0))*scaling*myλ0
            F!(F⁺,(myq1+q0)/2,(myq1-q0)/h,h/2)
            R[   1:nq]   .-= (h^2/2).*F⁺
            R[nq+1:nq+nλ] .= scaling*Φ(myq1)
        end
    end

    R_step1! = R_step1_maker(q0,q̇0,M,Φ,A,F!,nq,nλ,dt)
    initial_x .= vcat(q0,λ0)
    # @show q0
    # @show Φ(q0)
    # initial_R = similar(initial_x)
    # R_step1!(initial_R,initial_x)
    # @code_warntype R_step1!(initial_R,initial_x)
    # @show initial_R
    R_step1_result = nlsolve(R_step1!, initial_x, ftol=ftol)
    iteration = R_step1_result.iterations
    x1 = R_step1_result.zero
    push!(qs,x1[   1:nq])
    push!(q̇s,(qs[end] - q0)/dt)
    λs[1] .= x1[nq+1:nq+nλ]
    #---------Step 1 finisher-----------
    step += 1
    push!(ts,ts[end] + dt)
    sim.t = sim.tprev + dt
    sim.q .= qs[end]
    sim.q̇ .= q̇s[end]
    #---------Step 1 finisher-----------
    #---------Step 1 Callback-----------
    if callback.condition(sim)
        callback.affect!(sim)
    end
    #---------Step 1 Callback-----------

    function R_stepk_maker(qᵏ⁻¹,qᵏ,M,Φ,A,F!,nq,nλ,tᵏ,dt)
        # F!(F⁻,(qᵏ⁻¹.+qᵏ)./2,(qᵏ.-qᵏ⁻¹)./dt,tᵏ-dt/2)
        # kr = norm(F⁻,Inf)
        # scaling = mr + kr*dt^2
        function inner_R_stepk!(R,x)
            h = dt
            qᵏ⁺¹ = @view x[   1:nq]
            λᵏ   = @view x[nq+1:nq+nλ]
            R[   1:nq]    .= M*(qᵏ⁺¹.-2 .*qᵏ.+qᵏ⁻¹) .- transpose(A(qᵏ))*scaling*λᵏ
            F!(F⁺,(qᵏ⁺¹.+qᵏ)./2,(qᵏ⁺¹.-qᵏ)./h,tᵏ+h/2)
            F!(F⁻,(qᵏ⁻¹.+qᵏ)./2,(qᵏ.-qᵏ⁻¹)./h,tᵏ-h/2)
            R[   1:nq]   .-= h^2*1/2 .*( F⁺ .+ F⁻ )
            R[nq+1:nq+nλ] .= scaling*Φ(qᵏ⁺¹)
        end
    end

    if verbose
        @show iteration
    end
    #return SPARKState(ts,qs,q̇s,ps,λs)

    totalstep = ceil(Int,tspan[end]/dt)
    @showprogress 1 for timestep = 2:totalstep
        qᵏ⁺¹ = copy(qs[1])
        q̇ᵏ⁺¹ = copy(q̇s[1])
        λᵏ   = copy(λs[1])
        qᵏ⁻¹ = qs[end-1]
        qᵏ   = qs[end]
        λᵏ⁻¹ = λs[end]
        tᵏ   = ts[end]
        R_stepk! = R_stepk_maker(qᵏ⁻¹,qᵏ,M,Φ,A,F!,nq,nλ,tᵏ,dt)
        initial_x[   1:nq]    = qᵏ
        initial_x[nq+1:nq+nλ] = λᵏ⁻¹
        # initial_R = similar(initial_x)
        #R_stepk!(initial_R,initial_x)
        # @show initial_R
        #@code_warntype R_stepk!(initial_R,initial_x)
        R_stepk_result = nlsolve(R_stepk!, initial_x, ftol=ftol)
        #@code_warntype nlsolve(R_stepk!, initial_x, ftol=ftol)
        iteration += R_stepk_result.iterations
        xᵏ⁺¹ = R_stepk_result.zero
        qᵏ⁺¹ .= xᵏ⁺¹[   1:nq]
        λᵏ   .= xᵏ⁺¹[nq+1:nq+nλ]

        push!(qs,qᵏ⁺¹)
        push!(q̇s,(qᵏ⁺¹ - qᵏ)/dt)
        push!(λs,λᵏ)

        #---------Step k finisher-----------
        step += 1
        push!(ts,ts[end] + dt)
        sim.tprev = sim.t
        sim.qprev .= sim.q
        sim.t = sim.tprev + dt
        sim.q .= qs[end]
        sim.q̇ .= q̇s[end]
        #---------Step k finisher-----------
        #---------Step 1 Callback-----------
        if callback.condition(sim)
            callback.affect!(sim)
        end
        #---------Step 1 Callback-----------
        if verbose
            @printf("Progress: %5.1f%%, step: %s, time: %s, iterations: %s \n", (
            ts[end]/tspan[end]*100), timestep, ts[end], R_stepk_result.iterations)
        end
    end

    return sol
end
