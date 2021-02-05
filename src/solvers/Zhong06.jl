struct Zhong06Cache end
generate_cache(::Zhong06,intor) = Zhong06Cache()

function solve(intor::Integrator,cache::Zhong06Cache;dt=0.01,ftol=1e-14,verbose=false,
                                 callback=DEFAULT_CALLBACK)
    @unpack prob,sol = intor
    @unpack ts,qs,q̇s,ps,λs = sol
    @unpack funcs,tspan,q0,q̇0,λ0,nx,nq,nλ = prob
    M,Φ,A,F!,Jac_F! = funcs
    invM = inv(M)
    F⁺ = zero(q0)
    F⁻ = zero(q0)
    tend = tspan[2]
    step = 0
    initial_x = zeros(nx)
    initial_F = similar(initial_x)
    mr = norm(M,Inf)
    scaling = mr

    # function Momentum_k(qᵏ⁻¹,qᵏ,λᵏ,M,A,F!,tᵏ⁻¹,h)
    #     pᵏ = 1/h.*M*(qᵏ.-qᵏ⁻¹) .+ 1/h.*scaling.*transpose(A(qᵏ))*λᵏ
    #     F!(F⁺,(qᵏ.+qᵏ⁻¹)./2,(qᵏ.-qᵏ⁻¹)./h,tᵏ⁻¹+h/2)
    #     pᵏ .+= h/2 .*F⁺
    #     pᵏ
    # end
    function Momentum_k(qᵏ⁻¹,pᵏ⁻¹,qᵏ,λᵏ,M,A,h)
        pᵏ = -pᵏ⁻¹.+2/h.*M*(qᵏ.-qᵏ⁻¹) .+ 1/h.*scaling.*(
                transpose(A(qᵏ)) - transpose(A(qᵏ⁻¹)))*λᵏ
    end
    function R_stepk_maker(qᵏ⁻¹,pᵏ⁻¹,M,Φ,A,F!,nq,nλ,tᵏ⁻¹,dt)
        # F!(F⁻,(qᵏ⁻¹.+qᵏ)./2,(qᵏ.-qᵏ⁻¹)./dt,tᵏ-dt/2)
        # kr = norm(F⁻,Inf)
        # scaling = mr + kr*dt^2
        function inner_R_stepk!(R,x)
            h = dt
            qᵏ = @view x[   1:nq]
            λᵏ = @view x[nq+1:nq+nλ]
            R[   1:nq]    .= M*(qᵏ.-qᵏ⁻¹) .-h.*pᵏ⁻¹.- transpose(A(qᵏ⁻¹))*scaling*λᵏ
            F!(F⁺,(qᵏ.+qᵏ⁻¹)./2,(qᵏ.-qᵏ⁻¹)./h,tᵏ⁻¹+h/2)
            R[   1:nq]   .-= (h^2)/2 .*F⁺
            R[nq+1:nq+nλ] .= scaling*Φ(qᵏ)
        end
    end

    ∂F∂q = zeros(eltype(q0),nq,nq)
    ∂F∂q̇ = zeros(eltype(q0),nq,nq)

    function J_stepk_maker(qᵏ⁻¹,M,A,Jac_F!,nq,nλ,tᵏ⁻¹,dt)
        function inner_J_stepk!(J,x)
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
    totalstep = ceil(Int,tspan[end]/dt)
    @showprogress 1 for timestep = 1:totalstep
        #---------Step k Callback-----------
        if callback.condition(intor)
            callback.affect!(intor)
        end
        #---------Step k Callback-----------
        qᵏ = copy(qs[1])
        q̇ᵏ = copy(q̇s[1])
        pᵏ = copy(ps[1])
        λᵏ = copy(λs[1])
        qᵏ⁻¹ = qs[end]
        pᵏ⁻¹ = ps[end]
        λᵏ⁻¹ = λs[end]
        tᵏ⁻¹ = ts[end]
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
        R_stepk_result = nlsolve(dfk, initial_x, ftol=ftol, method = :newton)
        #@code_warntype nlsolve(R_stepk!, initial_x, ftol=ftol)
        if converged(R_stepk_result) == false
            error("Not Converged!")
        end
        iteration += R_stepk_result.iterations
        xᵏ = R_stepk_result.zero
        qᵏ .= xᵏ[   1:nq]
        λᵏ .= xᵏ[nq+1:nq+nλ]

        # pᵏ = Momentum_k(qᵏ⁻¹,qᵏ,λᵏ,M,A,F!,tᵏ⁻¹,dt)
        pᵏ = Momentum_k(qᵏ⁻¹,pᵏ⁻¹,qᵏ,λᵏ,M,A,dt)
        push!(qs,qᵏ)
        push!(ps,pᵏ)
        Aq = A(qᵏ)
        invM = inv(M+transpose(Aq)*Aq)
        push!(q̇s,invM*pᵏ)
        push!(λs,λᵏ)

        #---------Step k finisher-----------
        step += 1
        push!(ts,ts[end] + dt)
        intor.tprev = intor.t
        intor.qprev .= intor.q
        intor.q̇prev .= intor.q̇
        intor.t = intor.tprev + dt
        intor.q .= qs[end]
        intor.q̇ .= q̇s[end]
        #---------Step k finisher-----------
        if verbose
            @printf("Progress: %5.1f%%, step: %s, time: %s, iterations: %s \n", (
            ts[end]/tspan[end]*100), timestep, ts[end], R_stepk_result.iterations)
        end
    end

    return sol
end
