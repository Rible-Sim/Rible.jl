rO = 0.0025
rM = 0.0031
hM = 0.0058
lM = 0.010
lG = 0.015
hS = 0.02
lS = 0.0201
mM = 0.0003
mS = 0.0045
JM = 5.0*10^-9
JS = 7.0*10^-7

cᵩ = 0.0056
g = 9.81

εN = [0.5,0,0]
εT = [0,0,0]
μ = [0.3,0.3,0.3]

y0 = 0
φM0 = -0.1036
φS0 = -0.2788
ν0 = -0.3411
ωM0 = 0
ωS0 = -7.4583

q0 = [y0,φM0,φS0]
u0 = [ν0,ωM0,ωS0]
M = Matrix{Float64}(undef,3,3)
M[1,1] = mS + mM
M[1,2] = mS*lM
M[1,3] = mS*lG
M[2,1] = mS*lM
M[2,2] = JM + mS*lM^2
M[2,3] = mS*lM*lG
M[3,1] = mS*lG
M[3,2] = mS*lM*lG
M[3,3] = JS + mS*lG^2
M
h = Vector{Float64}(undef,3)
function hfunc!(h,q)
    y,φM,φS = q
    h[1] = -(mS+mM)*g
    h[2] = -cᵩ*(φM-φS) - mS*lM*g
    h[3] = -cᵩ*(φS-φM) - mS*lG*g
end
wN = [[0,0,-hS],
      [0, hM,0],
      [0,-hM,0]]
wT = [[1,lM,lG-lS],
      [1,rM,0],
      [1,rM,0]]

qA = q0
uA = u0
#qA = qs[it]
#uA = us[it]
dt = 0.01
qM = qA + 1/2*dt*uA
hfunc!(h,qM)
function gN1(q)
    y,φM,φS = q
    (lM + lG - lS - rO) - hS*φS
end
function gN2(q)
    y,φM,φS = q
    (rM - rO) + hM*φM
end
function gN3(q)
    y,φM,φS = q
    (rM - rO) - hM*φM
end
function indexset(q)
    H = Vector{Int64}()
    if gN1(q) <= 0
        push!(H,1)
    end
    if gN2(q) <= 0
        push!(H,2)
    end
    if gN3(q) <= 0
        push!(H,3)
    end
    H
end
𝓗 = indexset(qM)
indexset(q0)

M⁻¹ = inv(M)
hfunc!(h,qM)
WN = zeros(Float64,3,length(𝓗))
WT = zeros(Float64,3,length(𝓗))
μk = zeros(Float64,length(𝓗))
εNk = zeros(Float64,length(𝓗))
εTk = zeros(Float64,length(𝓗))
for (k,i) in enumerate(𝓗)
    WN[:,k] .= wN[i]
    WT[:,k] .= wT[i]
    μk[:,k] .= μ[i]
    εNk[:,k] .= εN[i]
    εNk[:,k] .= εN[i]
end
γANk = transpose(WN)*uA
γATk = transpose(WT)*uA
#uE = M⁻¹*(WN-WT*μ)*ΛN + M⁻¹*WT*ΛR + M⁻¹*h*dt + uA
E = Diagonal(ones(Float64,length(𝓗)))
A = [transpose(WN)*M⁻¹*(WN-WT*μk) transpose(WN)*M⁻¹*WT 0;
    transpose(WT)*M⁻¹*(WN-WT*μk) transpose(WT)*M⁻¹*WT E;
    2*μk                         -E                   0]
b = Vector{Float64}(undef,3)
b[1:1] .= transpose(WN)*M⁻¹*h*dt + (E + εNk)*γANk
b[2:2] .= transpose(WT)*M⁻¹*h*dt + (E + εTk)*γATk
b[3:3] .= 0

using PATHSolver



myfunc(x) = A*x + b

n = 3
lb = zeros(n)
ub = Inf*ones(n)

options(convergence_tolerance=1e-2, output=:yes, time_limit=3600)


z, f = solveLCP(myfunc, lb, ub)
