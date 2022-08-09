

function build_double()
    m1 = 2.6934977798E-02
    m2 = 0.1271425597
    l1 = 0.07071067811865477
    l2 = 2*0.052704627669472995
    I1 = 1.1222907415833337e-5
    I2 = 7.063475538888889e-5
    g = 9.81
    x1(θ1) =  l1/2*sin(θ1)
    y1(θ1) = -l1/2*cos(θ1)
    ẋ1(θ1,ω1) = l1/2*cos(θ1)*ω1
    ẏ1(θ1,ω1) = l1/2*sin(θ1)*ω1
    x2(θ1,θ2) =  l1*sin(θ1) + l2/2*sin(θ2)
    y2(θ1,θ2) = -l1*cos(θ1) - l2/2*cos(θ2)
    ẋ2(θ1,θ2,ω1,ω2) = l1*cos(θ1)*ω1 + l2/2*cos(θ2)*ω2
    ẏ2(θ1,θ2,ω1,ω2) = l1*sin(θ1)*ω1 + l2/2*sin(θ2)*ω2

    T(θ1,θ2,ω1,ω2) = (1//8*m1*l1^2 + 1//2*m2*l1^2 + 1//2*I1)*ω1^2 +
                     (1//8*m2*l2^2 + 1//2*I2               )*ω2^2 +
                      1//2*m2*cos(θ1-θ2)*l1*l2*ω1*ω2
    V(θ1,θ2) = m1*g*y1(θ1) + m2*g*y2(θ1,θ2)
    M(θ1,θ2) = [
        1//4*m1*l1^2 + m2*l1^2 + I1 1//2*m2*cos(θ1-θ2)*l1*l2;
        1//2*m2*cos(θ1-θ2)*l1*l2    1//4*m2*l2^2 + I2
    ]
    c(θ1,θ2) = 1//2*m2*sin(θ1-θ2)*l1*l2
    # function inner_f(ω,θ,p,t)
    #     ω1,ω2 = ω
    #     θ1,θ2 = θ
    #     -M(θ1,θ2)\[
    #          c(θ1,θ2)*ω2^2 + m1*g*l1/2*sin(θ1) + m2*g*l1*sin(θ1)
    #         -c(θ1,θ2)*ω1^2 + m2*g*l2/2*sin(θ2)
    #     ]
    # end
      function inner_f!(du,u,p,t)
        θ1,θ2 =  u[3:4]
        ω1,ω2 = du[3:4] .= M(θ1,θ2)\u[1:2]
        du[1:2] .= [
            -c(θ1,θ2)*ω1*ω2 - m1*g*l1/2*sin(θ1) - m2*g*l1*sin(θ1),
             c(θ1,θ2)*ω1*ω2 - m2*g*l2/2*sin(θ2)
        ]
    end
    inner_f!,M
end

function get_ref_sol(;saveat)
    θ_0 = [
        deg2rad(45), deg2rad(90)+0.32175055439664213
    ]
    ω_0 = zeros(2)
    u_0 = vcat(ω_0,θ_0)
    f!,M = build_double()
    prob = DE.ODEProblem(f!,u_0,(saveat[begin],saveat[end]))
    sol = DE.solve(prob,DE.Tsit5();dtmax=1e-5,abstol=1e-14,saveat)
    sol,M
end

function approx_slack(;k,μ,l0,filename="")
    μ⁺ = [-1e-1,-1e-2,-1e-3,-1e-4,-1e-5,-1e-6,-1e-7,-1e-8,-1e-9,0,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1]
    x = μ .+ μ⁺ .- l0
    y = [ifelse(xi+l0-μ>0,k*(xi+l0-μ),0) for xi in x]
    if filename != ""
        xy = Tables.table(hcat(x,y);header=["x",filename])
        CSV.write(filename*".csv",xy)
    end
    x,y
end

function random_vibration(;h=1e-3,tend=5.0)
    tend = 5.0
    N = round(Int,tend/h)
    t = [k*h for k = 0:N-1]
    f = OffsetArray([k/(N*h) for k = 1:N],0:N-1)
    psd = Uniform(minimum(f),maximum(f))
    # lines(psd)
    p = OffsetArray([pdf(psd,fk) for fk in f],0:N-1)
    a = OffsetArray(rand(Uniform(0,2π),N÷2+1),0:N÷2)
    A = OffsetArray(zeros(ComplexF64,N+1),0:N)
    for k = 0:N÷2
        A[k] = sqrt(p[k]/(2*N*h))*exp(im*a[k])
        A[N-k] = conj(A[k])
        if k in [0,N÷2]
            A[N-k] = A[k] = real(A[k])
        end
    end
    # A[N÷2]
    # x = [
    #     sum([
    #         A[(k-1)]*exp(-im*2π*(k-1)*(n-1)/N)
    #         for k = 1:N
    #     ]) |> real
    #     for n = 1:N
    # ]
    ẍ = fft(A[0:N-1]) .|> real
    # lines(x)
    ẋ = cumul_integrate(t,ẍ)
    # lines(y)
    x = cumul_integrate(t,ẋ)
    ẍ,ẋ,x
end
