
# plot(sol,vars=[3])
# @code_warntype f!(copy(u_0),u_0,0,0)
#
# prob_ode = SecondOrderODEProblem(f,ω_0,θ_0,(0.0,4.0))
# result = solve(prob_ode,DPRKN6();dtmax=1e-3)
#
# plot(result,vars = [3])
#
# T,V = build_double()
# V(0.1,0.1)
# T(0.1,0.1,0.1,0.1)
#
#
#
# using ModelingToolkit
# using Latexify
# using SymbolicUtils
# SymbolicUtils.show_simplified[] = true
# @parameters  m[1:2] l[1:2] α
# @variables t::Real θ[1:2](t)
# @variables ω[1:2], θ_[1:2]
# D = Differential(t)
# x1 =  l[1]/2*sin(θ[1])
# y1 = -l[1]/2*cos(θ[1])
# x2 = x1 + l[2]/2*sin(θ[2])
# y2 = y1 - l[2]/2*cos(θ[2])
# ẋ1 = expand_derivatives(D(x1))
# ẏ1 = expand_derivatives(D(y1))
# ẋ1 = substitute(ẋ1,Dict([D(θ[1]) => ω[1], D(θ[2]) => ω[2], θ[1] => θ_[1]]))
# ẏ1 = substitute(ẏ1,Dict([D(θ[1]) => ω[1], D(θ[2]) => ω[2], θ[1] => θ_[1]]))
# # r1 = @rule ~x*~y + ~x*~z => ~x*(~y+~z)
# simplify(ẋ1^2 + ẏ1^2; )
#
# ẋ2 = expand_derivatives(D(x2))
# ẏ2 = expand_derivatives(D(y2))
# ẋ2^2 + ẏ2^2
#
# using SymbolicUtils
#
# SymbolicUtils.show_simplified[] = true
#
# @syms x::Real y::Real
#
#
# 2x^2 - y + x^2
#
#
# x1 = (l[1]*sin(D(θ[1]))*ω[1])^2
# y1 = (l[1]*cos(D(θ[1]))*ω[1])^2
# x1 + y1
# istree(x1)
# simplify_fractions(x1 + y1)
# r = @rule sinh(im * ~x) => sin(~x)
#
#
# r(sinh(im * y))
#
#
# simplify(cos(y)^2 + sinh(im*y)^2, RuleSet([r]))
