using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using RecursiveArrayTools
using BenchmarkTools
using LaTeXStrings
using CSV
using GLMakie
using Revise
import Rible as RB
using JuMP
using Ipopt
using ForwardDiff
# import PyPlot as plt
cd("examples/manipulator")
includet("man_define.jl")
includet("man_plotting.jl")
includet("man_compare.jl")
includet("../analysis.jl")
# includet("two_sevor.jl")
# includet("../plot_helpers.jl")
# set_pyplot()
k = 250.0; c = 100.0;
# secs = 0; mics = 0
dt = 0.01
# dt = 0.01 # Same dt used for PID AND Dynamics solver
num_dof = 2
num_dof ==2 ? r̄ = [-0.28,0] : r̄ = [0.6,0.2]#3
#4-[0.45,-0.18,0.75,-0.25]
manipulator = man_ndof(num_dof;k,c)

function free(bot)
    st = bot.st
    M = RB.build_massmatrix(st)
    Φ = RB.build_Φ(st)
    A = RB.build_A(st)
    Q̃ = RB.build_Q̃(st)
    function F!(F,q,q̇,t)
        # u = -0.0025032clamp(t,0,1)
        # RB.actuate!(bot,[u])
        RB.reset_forces!(st)
        RB.distribute_q_to_rbs!(st,q,q̇)
        RB.update_cables_apply_forces!(st)
        RB.apply_gravity!(st)
        RB.assemble_forces!(F,st)
        # @show isapprox(F,Q̃*RB.fvector(st))
    end
    Jac_Γ = RB.build_Jac_Γ(st)
    function Jac_F!(∂F∂q,∂F∂q̇,q,q̇,t)
        ∂Γ∂q,∂Γ∂q̇ = Jac_Γ(q,q̇)
        Q̃*∂Γ∂q,Q̃*∂Γ∂q̇
    end
    M,Φ,A,F!,nothing
end

function make_control_lyapunov(bot)
    @unpack st = bot
    Q̃ = RB.build_Q̃(st)
    G = RB.build_G!(st)
    C3 = st.rigidbodies[end].state.cache.Cp[3]
    T_num_T = RB.build_Ti(st,st.nbodies)
    P = C3*T_num_T
    del_ang = zeros(Int64,num_servo,1)
    tar_angle_prev = init_pul
    #
        #=
        C2 = st.rigidbodies[end-2].state.cache.Cp[3]
        T2 =  RB.build_Ti(st,st.nbodies-2)
        C = zeros(4,8)
        T = zeros(8,st.ncoords)
        C[1:2,1:4] = C2
        C[3:4,5:8] = C3
        T[1:4,:] = T2
        T[5:8,:] = T_num_T
        P = C*T
        =#
        #=
        st = manipulator.st
        M = RB.build_massmatrix(st)
        invM = inv(M)
        q,q̇ = RB.get_coords(st)
        A = RB.build_A(st)
        Ak = A(q)
        Ȧk = A(q̇)
        transAk = transpose(Ak)
        C3 = st.rigidbodies[end].state.cache.Cp[3]
        T_num_T = RB.build_Ti(st,st.nbodies)
        P = C3*T_num_T
        =#

    function control_lyapunov!(intor,cache)
        @unpack prob,state,nx,nq,nλ = intor
        @unpack tspan,dyfuncs,restart = prob
        @unpack t,q,q̇,tprev,qprev,q̇prev = state
        @unpack totaltime,totalstep,ts,qs,q̇s,ps,λs,invM = cache
        M,Φ,A,F!,Jac_F! = dyfuncs
        # u = RB.get_cables_restlen(bot)
        # sum(u)
        # @show bot === prob.bot

        #=------BLF------

            RB.reset_forces!(st)
            RB.distribute_q_to_rbs!(st,q,q̇)
            RB.update_cables_apply_forces!(st)
            l = RB.get_cables_len(st)
            l̇ = RB.get_cables_len_dot(st)
            u_last = RB.get_cables_restlen(st)
            k = RB.get_cables_stiffness(st)
            c = RB.get_cables_damping(st)
            L = RB.build_L(st)
            𝛹 = I
            𝛩 = I
            Ak = A(q)
            transAk = transpose(Ak)
            Ȧk =  A(q̇)
            msn = inv(Ak*invM*transAk)
            Msn = invM*transAk*msn
            Nsn = invM - Msn*Ak*invM
            ϵ = 5/180*π

            function construct_s(r1,r2,r3)#::Vector{Float64}
                π/3 - acos(dot(r3-r2,r2-r1)/norm(r3-r2)/norm(r2-r1))
            end

            Ts1 = [1 0 0 0; 0 1 0 0]
            Ts2 = [0 0 1 0; 0 0 0 1]
            Tsi = zeros(2*st.nbodies+2,st.ncoords)
            for i in 1:1:st.nbodies
               Ti = RB.build_Ti(st,i)
               if i==1
                   Tsi[1:2,:] = Ts1*Ti
                   Tsi[3:4,:] = Ts2*Ti
               else
                   Tsi[i*2+1:i*2+2,:] = Ts2*Ti
               end
            end

            sx = zeros(st.num_of_dof)
            ∇s1 = zeros(st.num_of_dof,st.ncoords*2)
            for i in 1:1:st.num_of_dof
                    S(s::Vector) = construct_s(Tsi[i*2-1:i*2,:]*s,Tsi[i*2+1:i*2+2,:]*s,Tsi[i*2+3:i*2+4,:]*s)
                    s = q
                    sx[i] = S(s)
                    g = s -> ForwardDiff.gradient(S, s);
                    ∇s1[i,1:st.ncoords] = g(s)
            end
            sx
            V1(vx::Vector) = transpose(P*vx-r̄)*𝛩*(P*vx-r̄)
            V2(dvx::Vector) = transpose(P*dvx)*(P*dvx)
            vx = q
            dvx = q̇
            ∇V1 = vx -> ForwardDiff.gradient(V1, vx);
            ∇V2 = dvx -> ForwardDiff.gradient(V2, dvx);
            ∇V = zeros(1,2*st.ncoords)
            ∇V[1,1:st.ncoords] = ∇V1(vx)
            ∇V[1,st.ncoords+1:end] = ∇V2(dvx)
            fx = zeros(2*st.ncoords,1)
            fx[1:st.ncoords,1] = q̇
            fx[st.ncoords+1:end,1] = Nsn*G - Msn*Ȧk*q̇
            gx = zeros(2*st.ncoords,st.ncables)
            gx[st.ncoords+1:end,:] = Nsn*Q̃*L
            x = zeros(2*st.ncoords,1)
            x[1:st.ncoords,1] = q
            x[st.ncoords+1:end,1] = q̇


            if 0<sx[1]<ϵ && sx[2]>0 #&& sx[3]>0 #&& sx[4]>0 && sx[5]>0 && sx[6]>0
                ax = (∇V-((ϵ^2-s1[1]^2)/s1[1]^2).*∇s1[1,:])*fx
                bx = (∇V-((ϵ^2-s1[1]^2)/s1[1]^2).*∇s1[1,:])*gx
            elseif 0<sx[2]<ϵ && sx[1]>0 #&& sx[3]>0
                ax = (∇V-((ϵ^2-s1[2]^2)/s1[2]^2).*∇s1[2,:])*fx
                bx = (∇V-((ϵ^2-s1[2]^2)/s1[2]^2).*∇s1[2,:])*gx
            #elseif 0<sx[3]<ϵ && sx[1]>0 && sx[2]>0
                #ax = (∇V-((ϵ^2-s1[3]^2)/s1[3]^2).*∇s1[3,:])*fx
                #bx = (∇V-((ϵ^2-s1[3]^2)/s1[3]^2).*∇s1[3,:])*gx
            elseif sx[2]>ϵ && sx[1]>ϵ #&& sx[3]>ϵ
                ax = ∇V*fx
                bx = ∇V*gx
            end

            function η(α,β,μ)
                η = (-(norm(α)+sqrt(norm(α)^2)+μ*(norm(β))^4)/(norm(β))^2)*norm(β)
            end

            kx = 0.1+0.1*norm(x)
            μ  = kx/(kx+norm(bx)^2)

            if bx==0
                result_γ = zeros()
            else
                result_γ  = η(ax,bx,μ)
            end
            @show result_γ
            res_u = RB.force_densities_to_restlen(st,result_γ)
            for (s,u) in zip(st.cables,res_u)
                s.state.restlen = u
            end
            =#
            #------BLF†----------

        #-----Laypunov-----
            RB.reset_forces!(st)
            RB.distribute_q_to_rbs!(st,q,q̇)
            RB.update_cables_apply_forces!(st)
            l = RB.get_cables_len(st)
             #    l̇ = RB.get_cables_len_dot(st)

                # u_last = RB.get_cables_restlen(st)
                # k = RB.get_cables_stiffness(st)
                # c = RB.get_cables_damping(st)
                # a = 0.1;b = 0.08
                # l_min = sqrt(4*a^2-2*b^2+2*a*b)*1.05
                # l_max = (2*a + 2*b)*0.9
                # D = k - (u_last.*k-c.*l̇)./l_max

            #控制律
                L = RB.build_L(st)

                𝛹 = I
                𝛩 = I*10
                # transĖ = transpose(P*q̇)
                Ak = A(q)
                transAk = transpose(Ak)
                Ȧk =  A(q̇)
                msn = inv(Ak*invM*transAk)
                Msn = invM*transAk*msn
                Nsn = P*(invM - Msn*Ak*invM)
                L_u = diagm(l)
                invL_u = inv(L_u)
                L_u_diag = diag(invL_u)
                # T_u = [1 0;-1 0;0 1;0 -1]
                T_u = zeros(2*num_dof,num_dof)
                for i in 1:num_dof
                    T_u[2*i-1:2*i,i] = [1;-1]
                end
                u_0 = 0.18
                K_u = (-k*u_0).*L_u_diag.+k
                L_T = k.*invL_u*T_u

                ΓΓ = Nsn*Q̃*L*L_T
                μμ = Nsn*Q̃*L*K_u-P*Msn*Ȧk*q̇ + 𝛹*P*q̇ + 𝛩*(P*q-r̄) + Nsn*G  #Γγ<=μ

            #
            #barrier
                # ϵ = 40/180*π

                # function construct_s(r1,r2,r3)#::Vector{Float64}
                #     π/3 - acos(dot(r3-r2,r2-r1)/norm(r3-r2)/norm(r2-r1))
                # end

                # Ts1 = [1 0 0 0; 0 1 0 0]
                # Ts2 = [0 0 1 0; 0 0 0 1]
                # Tsi = zeros(2*st.nbodies+2,st.ncoords)
                # for i in 1:1:st.nbodies
                # Ti = RB.build_Ti(st,i)
                # if i==1
                #     Tsi[1:2,:] = Ts1*Ti
                #     Tsi[3:4,:] = Ts2*Ti
                # else
                #     Tsi[i*2+1:i*2+2,:] = Ts2*Ti
                # end
                # end

                # sx = zeros(st.num_of_dof)
                # ∇s1 = zeros(st.num_of_dof,st.ncoords)

                # for i in 1:1:st.num_of_dof
                #         S(s::Vector) = construct_s(Tsi[i*2-1:i*2,:]*s,Tsi[i*2+1:i*2+2,:]*s,Tsi[i*2+3:i*2+4,:]*s)
                #         s = q
                #         sx[i] = S(s)
                #         g = s -> ForwardDiff.gradient(S, s);
                #         ∇s1[i,:] = g(s)
                # end

                # if 0<sx[1]<ϵ && sx[2]>0 #&& sx[3]>0 #&& sx[4]>0 && sx[5]>0 && sx[6]>0
                #     μ =transĖ*P*Msn*Ȧk*q̇ - transĖ*𝛹*P*q̇ - transĖ*𝛩*(P*q-r̄) - transĖ*Nsn*G .+ (((ϵ^2-sx[1]^2)/sx[1]^2).*transpose(∇s1[1,:])*q̇)

                # elseif 0<sx[2]<ϵ && sx[1]>0 #&& sx[3]>0
                #     μ =transĖ*P*Msn*Ȧk*q̇ - transĖ*𝛹*P*q̇ - transĖ*𝛩*(P*q-r̄) - transĖ*Nsn*G .+ (((ϵ^2-sx[2]^2)/sx[2]^2).*transpose(∇s1[2,:])*q̇)

                # #elseif 0<sx[3]<ϵ && sx[1]>0 #&& sx[2]>0
                #     #μ =transĖ*P*Msn*Ȧk*q̇ - transĖ*𝛹*P*q̇ - transĖ*𝛩*(P*q-r̄) - transĖ*Nsn*G .+ (((ϵ^2-sx[3]^2)/sx[3]^2).*transpose(∇s1[3,:])*q̇)

                # elseif sx[2]>ϵ && sx[1]>ϵ #&& sx[3]>ϵ
                #     μ =transĖ*P*Msn*Ȧk*q̇ - transĖ*𝛹*P*q̇ - transĖ*𝛩*(P*q-r̄) - transĖ*Nsn*G

                # end
                #Nsn = P*(invM - Msn*Ak*invM)
                # Γ = transĖ*Nsn*Q̃*L
                #μ =transĖ*P*Msn*Ȧk*q̇ - transĖ*𝛹*P*q̇ - transĖ*𝛩*(P*q-r̄) - transĖ*Nsn*G#Γγ<=μ

            #Optimization
                Opt_Δu = Model(Ipopt.Optimizer)
                set_silent(Opt_Δu)
                @variable(Opt_Δu, Δu[1:num_dof])
                @objective(Opt_Δu, Min, sum(Δu))
                #@constraints(Opt_γ, con0, γ.>=0)
                @constraints(Opt_Δu, con0,Δu.<=0.04)
                @constraints(Opt_Δu, con1,Δu.>=-0.04)
                @constraints(Opt_Δu, con2, L_T*Δu.<=K_u)
                @constraints(Opt_Δu, con3, ΓΓ*Δu.==μμ)
                #@constraints(Opt_γ, con2, γ.<=D)
                optimize!(Opt_Δu)
                # @show termination_status(Opt_γ)
                #primal_status(Opt_γ)
                #dual_status(Opt_γ)
                objective_value(Opt_Δu)
                result_Δu = JuMP.value.(Δu)
                if result_Δu[1] >= 0.04
                    result_Δu[1] = 0.04
                elseif result_Δu[1]<=-0.04
                    result_Δu[1] = -0.04
                end
                if result_Δu[2] >= 0.04
                    result_Δu[2] = 0.04
                elseif result_Δu[2]<=-0.04
                    result_Δu[2] = -0.04
                end

                res_u = u_0 .+ T_u*result_Δu
                for (s,u) in zip(st.cables,res_u)
                    s.state.restlen = u
                end
                E = P*q-r̄

                for i in 1:num_dof
                    del_ang[i] = floor(Int,180/pi/r_rot*2000/270*result_Δu[i])
                    # if iseven(i)
                    #     del_ang[i] = -del_ang[i]
                    # end
                end
                tar_angle = init_pul + del_ang
                if tar_angle[1]>2400
                    tar_angle[1]=2400
                elseif tar_angle[1]<600
                    tar_angle[1]=600
                end
                if tar_angle[2]>2400
                    tar_angle[2]=2400
                elseif tar_angle[2]<600
                    tar_angle[2]=600
                end
                diff_tar_angle = tar_angle - tar_angle_prev
                diff_tar_angle[1] = abs(diff_tar_angle[1])
                diff_tar_angle[2] = abs(diff_tar_angle[2])
                mov_time = findmax(diff_tar_angle)[1]
                if mov_time<=20
                    mov_time=20
                elseif mov_time>=30000
                    mov_time=30000
                end
                servo_run(tar_angle,mov_time)
                LP.time_sleep(mov_time/1000)
                # angle = for_lya_control(result_Δu)
                # q = angle2q(angle)
                starttick = gpioTick()
                endtick = gpioTick()
                difftick_mics = convert(Int64, endtick - starttick)
                difftick_sec = difftick_mics/1000000
                @show res_u#l̇

    end
end

function lya_sevor_run(dt,t_tol,manipulator)
    st = manipulator.st
    C3 = st.rigidbodies[end].state.cache.Cp[3]
    T_num_T = RB.build_Ti(st,st.nbodies)
    P = C3*T_num_T

    init_pul = [1500,1500]
    result_Δu = zeros(num_dof,1)
    rot_init = zeros(Int16, 1,num_dof)
    rot_data = zeros(Int16, 1,num_dof)
    del_ang = zeros(Int64,num_servo,1)
    T_u = zeros(2*num_dof,num_dof)
    q̇ = zeros(2*num_dof+4,1)
    i2cReadI2CBlockData(i2c_handle,0x00,rot_init,2*num_dof)
    q_angle_0 = sensor_angle(rot_init)
    q = angle2q(q_angle_0)
    q_prev = q
    tar_angle_prev = init_pul
    𝛹 = I*0.05
    𝛩 = I*0.05
    starttick = gpioTick()

    for step in 1:dt:t_tol
        RB.reset_forces!(st)
        RB.distribute_q_to_rbs!(st,q,q̇)
        RB.update_cables_apply_forces!(st)
        Q̃ = RB.build_Q̃(st)
        # G = RB.build_G!(st)
        A = RB.build_A(st)
        M = RB.build_massmatrix(st)
        invM = inv(M)
        l = RB.get_cables_len(st)
        L = RB.build_L(st)
        G = RB.build_G!(st)
        #控制律
            Ak = A(q)
            transAk = transpose(Ak)
            Ȧk =  A(q̇)
            msn = inv(Ak*invM*transAk)
            Msn = invM*transAk*msn
            Nsn = P*(invM - Msn*Ak*invM)
            L_u = diagm(l)
            invL_u = inv(L_u)
            L_u_diag = diag(invL_u)
                for i in 1:num_dof
                    T_u[2*i-1:2*i,i] = [1;-1]
                end
            u_0 = 0.18
            K_u = (-k*u_0).*L_u_diag.+k
            L_T = k.*invL_u*T_u

            Γ = Nsn*Q̃*L*L_T
            μ = Nsn*Q̃*L*K_u - P*Msn*Ȧk*q̇ + 𝛹*P*q̇ + 𝛩*(P*q-r̄) + Nsn*G
        #Optimization
            Opt_Δu = Model(Ipopt.Optimizer)
            set_silent(Opt_Δu)
            @variable(Opt_Δu, Δu[1:num_dof])
            @objective(Opt_Δu, Min, sum(Δu))
            @constraints(Opt_Δu, con0,Δu.<=0.04)
            @constraints(Opt_Δu, con1,Δu.>=-0.04)
            @constraints(Opt_Δu, con2, L_T*Δu.<=K_u)
            @constraints(Opt_Δu, con3, Γ*Δu.==μ)
            optimize!(Opt_Δu)
            objective_value(Opt_Δu)
            result_Δu = JuMP.value.(Δu)
            if result_Δu[1] >= 0.04
                result_Δu[1] = 0.04
            elseif result_Δu[1]<=-0.04
                result_Δu[1] = -0.04
            end
            if result_Δu[2] >= 0.04
                result_Δu[2] = 0.04
            elseif result_Δu[2]<=-0.04
                result_Δu[2] = -0.04
            end

        #Opt_result
            res_u = u_0 .+ T_u*result_Δu
            for (s,u) in zip(st.cables,res_u)
                s.state.restlen = u
            end
            E = P*q-r̄
            traj = P*q
        #control
            for i in 1:num_dof
                del_ang[i] = floor(Int,180/pi/r_rot*2000/270*result_Δu[i])
                # if iseven(i)
                #     del_ang[i] = -del_ang[i]
                # end
            end
            tar_angle = init_pul + del_ang
            if tar_angle[1]>2400
                tar_angle[1]=2400
            elseif tar_angle[1]<600
                tar_angle[1]=600
            end
            if tar_angle[2]>2400
                tar_angle[2]=2400
            elseif tar_angle[2]<600
                tar_angle[2]=600
            end
            diff_tar_angle = tar_angle - tar_angle_prev
            diff_tar_angle[1] = abs(diff_tar_angle[1])
            diff_tar_angle[2] = abs(diff_tar_angle[2])
            mov_time = 2*findmax(diff_tar_angle)[1]
            if mov_time<=20
                mov_time=20
            elseif mov_time>=30000
                mov_time=30000
            end
            servo_run(tar_angle,mov_time)
            LP.time_sleep(mov_time/1000)
            i2cReadI2CBlockData(i2c_handle,0x00,rot_data,2*num_dof)
            angle = sensor_angle(rot_data)
            q = angle2q(angle)
            # angle = for_lya_control(result_Δu)
            endtick = gpioTick()
            difftick_mics = convert(Int64, endtick - starttick)
            difftick_sec = difftick_mics/1000000
            q̇ = (q - q_prev)/difftick_sec

            q_prev = q
            starttick = endtick
            tar_angle_prev = tar_angle
            # @show difftick_sec
            # @show mov_time*1e-3
            @show angle
            # @show q
    end
end

t_tol =1.5
lya_sevor_run(dt,t_tol,manipulator)

prob = RB.SimProblem(manipulator,free,make_control_lyapunov(manipulator),(0.0,10.0))
RB.solve!(prob,RB.Zhong06();dt,ftol=1e-14,exception=false)

prob = RB.SimProblem(manipulator,free,(0.0,10.0))
RB.solve!(prob,RB.Zhong06();dt,ftol=1e-14,exception=false)

prob = RB.SimProblem(man_inv,free,(0.0,10.0))
RB.solve!(prob,RB.Zhong06();dt,ftol=1e-14,exception=false)

plotstructure(manipulator)
using RowEchelon

A = RB.build_A(manipulator.st)
Q̃ = RB.build_Q̃(manipulator.st)
L̂ = RB.build_L̂(manipulator.st)
q,_ = RB.get_coords(manipulator.st)
A(q)
rref(hcat(transpose(A(q)),-Q̃*L̂))

get_angles(manipulator)
# RB.get_cables_deform(manipulator)
f1,f2 = RB.get_cables_tension(manipulator)
f2/f1

servo_run(init_pul,2000) #复位
singal_servo_run(1,[600,600],4000)

# rot_data = zeros(Int16, 1,2)
# i2cReadI2CBlockData(i2c_handle,0x00,rot_data,4)
# rot_data
rot_data = zeros(Int16, 1,num_dof)
f = open("freewave2.txt","w")
for i in 1:10000
    i2cReadI2CBlockData(i2c_handle,0x00,rot_data,2*num_dof)
    angle = sensor_angle(rot_data)
    println(f,angle[1])

end
close(f)

for_lya_control([-0.03,0])
