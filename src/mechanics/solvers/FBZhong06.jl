

struct FBZhong06Cache{MMT, AT, ΦT, ΨT}
    mass_matrices::MMT
    A::AT
    Φ::ΦT
    Ψ::ΨT
end

function generate_cache(
        simulator::Simulator{DynamicsProblem{RobotType,Contactless,SlidingTendon}},
        solver::DynamicsSolver{Zhong06};
        dt,kargs...
    )   where RobotType
    (;prob) = simulator
    (;bot,dynfuncs) = prob
    # F!,_ = dynfuncs
    mm = build_mass_matrices(bot)
    # (;M) = mm
    A = make_cstr_jacobian(bot)
    Φ = make_cstr_function(bot)
    Ψ = make_Ψ(bot.structure)
    FBZhong06Cache(mm,A,Φ,Ψ)
end

function solve!(simulator::Simulator,cache::FBZhong06Cache;
                dt,ftol=1e-14,verbose=false,iterations=50,
                progress=true,exception=true)
    (;prob,controller,tspan,restart,totalstep) = simulator
    (;bot,dynfuncs) = prob
    (;traj) = bot
    (;F!,Jac_F!) = dynfuncs
    (;actuate!) = controller
    (;mass_matrices,A,Φ,Ψ) = cache
    (;Ḿ,M̌,M̄,M̌⁻¹) = mass_matrices
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
            ζ = build_ζ(bot.structure)
            ∂ζ∂q = build_∂ζ∂q(bot.structure, q̌ᵏ)
            ∂ζ∂s̄ = build_∂ζ∂s̄(bot.structure)
            n = length(ζ)
            ∂s̄∂s̄ = I(n)
            κ₁ = 10; κ₂ = 10
            coes = diagm([((ζ[i]/κ₁)^2 + (κ₂*sᵏ[i])^2)^(-1/2) for i in 1:n])
            coζ = coes*diagm([ζ[i]/κ₁^2 for i in 1:n]) - diagm([1/κ₁ for i in 1:n])
            cos̄ = coes*diagm([κ₂^2*sᵏ[i] for i in 1:n]) - diagm([κ₂ for i in 1:n])
            # Jac_F!(∂F∂q̌,∂F∂q̌̇,(qᵏ.+qᵏ⁻¹)./2,(qᵏ.-qᵏ⁻¹)./h,tᵏ⁻¹+h/2)
            ∂Q̌∂q̌̇ = build_∂Q̌∂q̌̇(bot.structure)
            ∂Q̌∂q̌ = build_∂Q̌∂q̌(bot.structure)
            ∂Q̌∂s̄ = build_∂Q̌∂s̄(bot.structure)
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
            ζ = build_ζ(bot.structure)
            ∂ζ∂q = Record_build_∂ζ∂q(bot.structure, q̌ᵏ,"t1.xlsx","t1")
            ∂ζ∂s̄ = build_∂ζ∂s̄(bot.structure)
            n = length(ζ)
            ∂s̄∂s̄ = I(n)
            κ₁ = 10; κ₂ = 10
            coes = diagm([((ζ[i]/κ₁)^2 + (κ₂*sᵏ[i])^2)^(-1/2) for i in 1:n])
            coζ = coes*diagm([ζ[i]/κ₁^2 for i in 1:n]) - diagm([1/κ₁ for i in 1:n])
            cos̄ = coes*diagm([κ₂^2*sᵏ[i] for i in 1:n]) - diagm([κ₂ for i in 1:n])
            # Jac_F!(∂F∂q̌,∂F∂q̌̇,(qᵏ.+qᵏ⁻¹)./2,(qᵏ.-qᵏ⁻¹)./h,tᵏ⁻¹+h/2)
            ∂Q̌∂q̌̇ = build_∂Q̌∂q̌̇(bot.structure)
            ∂Q̌∂q̌ = build_∂Q̌∂q̌(bot.structure)
            ∂Q̌∂s̄ = build_∂Q̌∂s̄(bot.structure)
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
        # control!(sim,cache)
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
        if !(actuate! isa Missing)
            actuate!(bot.structure, tᵏ⁻¹; dt)
        end
        Res_stepk! = make_Res_stepk(qᵏ,q̌ᵏ,λᵏ,sᵏ,qᵏ⁻¹,p̌ᵏ⁻¹,F̌,Aᵀ,tᵏ⁻¹)
        if Jac_F! isa Missing
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
                # sim.convergence = false
                break
            end
        end
        total_iterations += Res_stepk_result.iterations
        xᵏ = Res_stepk_result.zero
        q̌ᵏ .= xᵏ[   1:nq̌   ]
        λᵏ .= xᵏ[nq̌+1:nq̌+nλ]
        sᵏ .= xᵏ[nq̌+nλ+1:nx]
        Momentum_k!(p̌ᵏ,p̌ᵏ⁻¹,qᵏ,qᵏ⁻¹,λᵏ,Ḿ,A,Aᵀ,dt)
        q̌̇ᵏ .= M̌⁻¹*(p̌ᵏ.-M̄*q̃̇ᵏ )
        #---------Step k finisher-----------
        step += 1
        p̌ᵏ,p̌ᵏ⁻¹ = p̌ᵏ⁻¹,p̌ᵏ
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

    bot
end
