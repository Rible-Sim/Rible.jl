using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using TableView
# using BenchmarkTools
# import PyPlot; const plt = PyPlot
# using LaTeXStrings
# using NLsolve
using HomotopyContinuation
using DynamicPolynomials
using Revise
using TensegritySolvers; const TS = TensegritySolvers
using TensegrityRobot
const TR = TensegrityRobot
include("tail_define.jl")

n = 1
tail = make_tail(n)
# @code_warntype make_tail(n)
q0,q̇0,λ0 = TR.get_initial(tail)

function dynfuncs(tgstruct,q0)

    M = TR.build_massmatrix(tgstruct)
    Φ = TR.build_Φ(tgstruct,q0)
    A = TR.build_A(tgstruct)

    #Q̃=TR.build_Q̃(tgstruct)

    function F!(F,q,q̇,t)
        TR.reset_forces!(tgstruct)
        TR.distribute_q_to_rbs!(tgstruct,q,q̇)
        TR.update_strings_apply_forces!(tgstruct)
        # F .= Q̃*TR.fvector(tgstruct)
        F .= 0.0
        TR.assemble_forces!(F,tgstruct)
    end

    M,Φ,A,F!,nothing
end

M,Φ,A,F!,Jacs = dynfuncs(tail,q0)

@var q[1:tail.ncoords]
@var λ[1:tail.nconstraint]
@var u[1:tail.nstrings]
@var a[1:tail.nstrings]
polyq = (1.0q)
polyλ = (1.0λ)
polyu = (1.0u)
polya = (1.0a)

Q̃ = TR.build_Q̃(tail)

U = TR.build_U(tail,polya,polyu)


S = TR.build_S(tail,polya,polyq)
TR.update_strings_apply_forces!(tail)
l = [s.state.length for s = tail.strings]
a0 = 1 ./l
# S[1](a=>a0,q=>q0)
# S[2](a=>a0,q=>q0)
# S[3](a=>a0,q=>q0)
# S[4](a=>a0,q=>q0)

@var g
G = TR.build_G(tail)
Fu = [-transpose(A(polyq))*polyλ  + g*G + Q̃*U*polyq;
    S;
    Φ(polyq)]
F =  [subs(f, u=>l) for f in Fu]

[f(q=>q0,λ=>λ0,a=>1 ./l) for f in F]
startsols = [[q0;λ0; 1 ./l]]
vg = [q,λ,a]
Fsys = System(Fu; variable_groups=vg, parameters = [u;g])

u0 = l
u1 = l - [0.0,0.0,0.01,0.0]
g0 = 0.0
g1 = 1.0
startparas = vcat(u0,g0)
endparas = vcat(u1,g1)
# F =  [subs(f, q=>q0,λ=>λ0,a=>1 ./l,u=>u0) for f in Fu]
# paths_to_track(Fsys; start_system = :total_degree)
#
# paths_to_track(Fsys; start_system = :polyhedral)

result = solve(Fsys, startsols; start_parameters=startparas, target_parameters=endparas)
results(result;only_real = true)
qλa = real_solutions(result)
q0

Fexp = expressions(Fsys)
@time solve(Fexp,show_progress=false)

@var x y z p[1:3]

F = System(
    [
        x + 3 + 2y + 2y^2 - p[1],
        (x - 2 + 5y) * z + 4 - p[2] * z,
        (x + 2 + 4y) * z + 5 - p[3] * z,
    ];
    parameters = p
)

p₀ = randn(ComplexF64, 3)
# Compute all solutions for F_p₀
result_p₀ = solve(F, target_parameters = p₀)
data = [randn(3) for _ in 1:10_000]
data_points = solve(
    F,
    solutions(result_p₀);
    start_parameters =  p₀,
    target_parameters = data
)

data_points = solve(
    F,
    solutions(result_p₀);
    start_parameters =  p₀,
    target_parameters = data,
    transform_result = (r,p) -> real_solutions(r),
    flatten = true
)


@var z[1:6]
α = randn(5)
a = randn(9)
# define the system of polynomials
f = [z[i,:] ⋅ z[i,:] for i = 2:5]
Expression.(z)
@polyvar z[1:6,1:3]
vg=[[z[2,:]; z[4,:]], [z[3,:]; z[5,:]]]
myvg =
@polyvar x y # assigns x (resp. y) to a variable of name x (resp. y)
p = 2x + 3.0x*y^2 + y
@test differentiate(p, x) # compute the derivative of p with respect to x
@test differentiate.(p, (x, y)) # compute the gradient of p
@test p((x, y)=>(y, x)) # replace any x by y and y by x
@test subs(p, y=>x^2) # replace any occurence of y by x^2
@test p(x=>1, y=>2) # evaluate p at [1, 2]

n = 3
A = rand(n, n)
@polyvar x[1:n] # assign x to a tuple of variables x1, x2, x3
p = sum(x .* x) # x_1^2 + x_2^2 + x_3^2
subs(p, x[1]=>2, x[3]=>3) # x_2^2 + 13
p(x=>A*vec(x))

Φrr = tail.rigidbodies[1].state.cache.funcs.Φ
Φqrr = tail.rigidbodies[1].state.cache.funcs.Φq
@polyvar q[1:4]
@polyvar λ[1:1]
Φqrr(polynomial.(1.0q))
ArrT = transpose(Φqrr(polynomial.(1.0q)))
ArrT*λ
PolyΦrr = Φrr(polynomial.(1.0q))
differentiate.(PolyΦrr, q)
Matrix{eltype(q)}(undef,1,2)
polynomial(q)
