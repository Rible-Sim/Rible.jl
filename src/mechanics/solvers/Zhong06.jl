struct Zhong06 <: AbstractSolver end

struct Zhong06Cache{solT,MMT,AT,ΦT}
    solver::solT
    mass_matrices::MMT
    A::AT
    Φ::ΦT
end

function generate_cache(solver::Zhong06,intor;dt,kargs...)
    (;prob) = intor
    (;bot,dynfuncs) = prob
    # F!,_ = dynfuncs
    mm = build_mass_matrices(bot)
    # (;M) = mm
    A = make_cstr_jacobian(bot)
    Φ = make_cstr_function(bot)
    Zhong06Cache(solver,mm,A,Φ)
end

function retrieve!(intor,cache::Zhong06Cache)

end

function solve!(intor::Integrator,cache::Zhong06Cache;
                dt,ftol=1e-14,verbose=false,maxiters=50,
                progress=true,exception=true)
    (;prob,controller,tspan,restart,totalstep) = intor
    (;bot,dynfuncs) = prob
    (;traj) = bot
    (;F!,Jac_F!) = dynfuncs
    (;mass_matrices,A,Φ) = cache
    (;Ḿ,M̌,M̄,invM̌) = mass_matrices
    q̌0 = traj.q̌[begin]
    λ0 = traj.λ[begin]
    q̇0 = traj.q̇[begin]
    T = eltype(q̌0)
    nq̌ = length(q̌0)
    nλ = length(λ0)
    ∂F∂q̌ = zeros(T,nq̌,nq̌)
    ∂F∂q̌̇ = zeros(T,nq̌,nq̌)
    p̌ᵏ⁻¹ = Ḿ*q̇0
    p̌ᵏ   = zero(p̌ᵏ⁻¹)
    step = 0
    initial_x = vcat(q̌0,λ0)
    initial_Res = zero(initial_x)
    initial_Jac = initial_Res*transpose(initial_Res)
    mr = norm(Ḿ,Inf)
    scaling = mr

    function make_Res_stepk(qᵏ,q̌ᵏ,λᵏ,qᵏ⁻¹,p̌ᵏ⁻¹,F̌,Aᵀ,tᵏ⁻¹)
        @inline @inbounds function inner_Res_stepk!(Res,x)
            h = dt
            q̌ᵏ .= x[   1:nq̌   ]
            λᵏ .= x[nq̌+1:nq̌+nλ]
            F!(F̌,(qᵏ.+qᵏ⁻¹)./2,(qᵏ.-qᵏ⁻¹)./h,tᵏ⁻¹+h/2)
            Res[   1:nq̌   ] .= Ḿ*(qᵏ.-qᵏ⁻¹) .-
                               h.*p̌ᵏ⁻¹ .-
                               (h^2)/2 .*F̌ .+
                               scaling.*Aᵀ*λᵏ
            Res[nq̌+1:nq̌+nλ] .= scaling.*Φ(qᵏ)
        end
    end

    function make_Jac_stepk(qᵏ,q̌ᵏ,qᵏ⁻¹,∂F∂q̌,∂F∂q̌̇,Aᵀ,tᵏ⁻¹)
        @inline @inbounds function inner_Jac_stepk!(Jac,x)
            h = dt
            q̌ᵏ .= x[1:nq̌]
            Jac_F!(∂F∂q̌,∂F∂q̌̇,(qᵏ.+qᵏ⁻¹)./2,(qᵏ.-qᵏ⁻¹)./h,tᵏ⁻¹+h/2)
            Jac[   1:nq̌ ,   1:nq̌ ] .=  M̌.-(h^2)/2 .*(1/2 .*∂F∂q̌.+1/h.*∂F∂q̌̇)
            Jac[   1:nq̌ ,nq̌+1:end] .=  scaling.*Aᵀ
            Jac[nq̌+1:end,   1:nq̌ ] .=  scaling.*A(qᵏ)
            Jac[nq̌+1:end,nq̌+1:end] .=  0.0
        end
    end

    @inline @inbounds function Momentum_k!(p̌ᵏ,p̌ᵏ⁻¹,qᵏ,qᵏ⁻¹,λᵏ,Ḿ,A,Aᵀ,h)
        p̌ᵏ .= -p̌ᵏ⁻¹.+2/h.*Ḿ*(qᵏ.-qᵏ⁻¹) .-
            1/h.*scaling.*(transpose(A(qᵏ))-Aᵀ)*λᵏ
    end

    iteration_break = 0
    prog = Progress(totalstep; dt=1.0, enabled=progress)
    dg_step = ceil(Int,log10(totalstep))+1
    dg_dt = max(1,-floor(Int,log10(dt)))
    wd_t = ceil(Int,log10(traj.t[end]))+dg_dt+1+1
    progfmt = Printf.Format("Progress: %5.1f%%, step: %$(dg_step)u, time: %$(wd_t).$(dg_dt)f, iterations: %s \n")
    
    for timestep = 1:totalstep
        #---------Step k Control-----------
        # control!(intor,cache)
        #---------Step k Control-----------
        tᵏ⁻¹ = traj.t[timestep]
        tᵏ = traj.t[timestep+1]
        qᵏ⁻¹ = traj.q[timestep]
        q̌ᵏ⁻¹ = traj.q̌[timestep]
        q̌̇ᵏ⁻¹ = traj.q̌̇[timestep]
        qᵏ = traj.q[timestep+1]
        q̌ᵏ = traj.q̌[timestep+1]
        q̌̇ᵏ = traj.q̌̇[timestep+1]
        q̃̇ᵏ = traj.q̃̇[timestep+1]
        λᵏ = traj.λ[timestep+1]
        F̌ = traj.F̌[timestep+1]
        initial_x[   1:nq̌]    .= q̌ᵏ⁻¹ .+ dt.*q̌̇ᵏ⁻¹
        initial_x[nq̌+1:nq̌+nλ] .= 0
        Aᵀ = transpose(A(qᵏ⁻¹))
        Res_stepk! = make_Res_stepk(qᵏ,q̌ᵏ,λᵏ,qᵏ⁻¹,p̌ᵏ⁻¹,F̌,Aᵀ,tᵏ⁻¹)
        isconverged = false
        if Jac_F! isa Nothing
            dfk = OnceDifferentiable(Res_stepk!,initial_x,initial_Res)
            Res_stepk_result = nlsolve(dfk, initial_x; ftol, iterations=maxiters, method=:newton)
            isconverged = converged(Res_stepk_result)
            iteration_break = Res_stepk_result.iterations
            Res_stepk_result.iterations
            xᵏ = Res_stepk_result.zero
        else
            Jac_stepk! = make_Jac_stepk(qᵏ,q̌ᵏ,qᵏ⁻¹,∂F∂q̌,∂F∂q̌̇,Aᵀ,tᵏ⁻¹)
            # if timestep == 10
            #     Jac_ref = zeros(nq̌+nλ,nq̌+nλ)
            #     FiniteDiff.finite_difference_jacobian!(Jac_ref,Res_stepk!,initial_x)
            #     Jac_my = zeros(nq̌+nλ,nq̌+nλ)
            #     Jac_stepk!(Jac_my,initial_x)
            #     diff_Jac = Jac_my .- Jac_ref
            #     diff_Jac[abs.(diff_Jac).<1e-5] .= 0.0
            #     diff_rowidx = findall((x)-> x>1e-5,norm.(eachrow(diff_Jac)))
            #     @show diff_rowidx .- length(q̌ᵏ)                
            #     @show diff_Jac[end-4:end-3,:]
            #     @show findall((x)->abs(x)>0,diff_Jac[114+length(q̌ᵏ),:])
            #     @show findall((x)->abs(x)>0,diff_Jac[115+length(q̌ᵏ),:])
            #     @show Jac_my[end-4:end-3,:]
            #     @show Jac_ref[end-4:end-3,:]
            # end
            # @show maximum(abs.(diff_Jac))for iteration = 1:maxiters
            # @show iteration,D,ηs,restitution_coefficients,gaps
            for iteration = 1:maxiters
                Res_stepk!(initial_Res,initial_x)
                normRes = norm(initial_Res)
                if normRes < ftol
                    isconverged = true
                    iteration_break = iteration-1
                    break
                end                
                # FiniteDiff.finite_difference_jacobian!(initial_Jac,Res_stepk!,initial_x,Val{:central})
                Jac_stepk!(initial_Jac,initial_x)
                initial_x .+= -initial_Jac\initial_Res
            end
            xᵏ = initial_x
            # dfk = OnceDifferentiable(Res_stepk!,Jac_stepk!,initial_x,initial_Res)
        end
        
        if !isconverged
            if exception
                error("Not Converged! Step=$timestep")
            else
                # intor.convergence = false
                break
            end
        end
        q̌ᵏ .= xᵏ[   1:nq̌   ]
        λᵏ .= xᵏ[nq̌+1:nq̌+nλ]
        Momentum_k!(p̌ᵏ,p̌ᵏ⁻¹,qᵏ,qᵏ⁻¹,λᵏ,Ḿ,A,Aᵀ,dt)
        q̌̇ᵏ .= invM̌*(p̌ᵏ.-M̄*q̃̇ᵏ )
        #---------Step k finisher-----------
        step += 1
        p̌ᵏ,p̌ᵏ⁻¹ = p̌ᵏ⁻¹,p̌ᵏ
        #---------Step k finisher-----------
        if verbose
            progstr = Printf.format(
                progfmt,
                floor(timestep/totalstep*100;digits=1), timestep, tᵏ, iteration_break
            )
            print(progstr)
        end
        next!(prog)
    end

    bot
end
