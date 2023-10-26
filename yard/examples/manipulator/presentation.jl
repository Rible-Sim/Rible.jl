using LinearAlgebra
using SparseArrays
using Parameters
using StaticArrays
using Makie
import PyPlot as plt
plt.pygui(true)
using LaTeXStrings
using RecursiveArrayTools
using Interpolations
using BenchmarkTools
using Revise
import Rible as RB

cd("examples/manipulator")
includet("man_define.jl")
includet("man_plotting.jl")
includet("man_compare.jl")
includet("../analysis.jl")
includet("../plot_helpers.jl")
set_pyplot()
includet("IScontrol.jl")

k = 1.e2
c = 50.0
target_t = 20.0
tend = 40.0

man_linear, aitp_linear = simulate_linearactuating(;k,c,target_t,tend,Alg=Linear())
man_quadra, aitp_quadra = simulate_linearactuating(;k,c,target_t,tend,Alg=Quadratic(Flat(OnGrid())))

tstops = [0,5,10,15,20,25]
fig = pyplot2(man_linear,tstops)

fig = plot_restlen(man_linear,aitp_linear,aitp_quadra)

plotstructure(man_linear)
