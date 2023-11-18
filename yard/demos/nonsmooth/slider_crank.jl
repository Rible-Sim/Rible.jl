figdir::String = ""
if Sys.iswindows()
    figdir::String = raw"C:\Users\luo22\OneDrive\Papers\FrictionalContact\CMAME"
elseif Sys.isapple()
    figdir::String = raw"/Users/jacob/Library/CloudStorage/OneDrive-SharedLibraries-onedrive/Papers/FrictionalContact/CMAME"
end
#-- point mass end

#--  Spinning top
include("deps.jl")
using AbbreviatedStackTraces #jl
ENV["JULIA_STACKTRACE_ABBREVIATED"] = true #jl
ENV["JULIA_STACKTRACE_MINIMAL"] = true #jl
import Rible as RB
include("../../vis.jl")
includet("../../vis.jl") #jl
include("../../../examples/robots/slider_crank.jl")
includet("../../../examples/robots/slider_crank.jl")#jl
sc = slider_crank(;coordsType=RB.QCF.QC)
plot_traj!(sc;showmesh=false,showground=false)
dt = 1e-3
tspan = (0.0,1.0)

function sc_contact_dynfuncs(bot;)
    ## horizontal Plane
    planes = [
        RB.Plane([0,0, 1.0],[0,0,-0.026]),
        RB.Plane([0,0,-1.0],[0,0, 0.026])
    ]
    RB.frictionless_contact_dynfuncs(bot;flatplane = planes)
end

prob = RB.SimProblem(sc,(bot)->RB.dynfuncs(bot;gravity=true))
RB.solve!(
    prob,
    RB.Zhong06();
    dt,tspan,ftol=1e-14,maxiters=50,verbose=true,exception=true,progress=false,
)

prob = RB.SimProblem(sc,sc_contact_dynfuncs)
RB.solve!(
    prob,
    RB.ZhongQCCPN();
    dt,tspan,ftol=1e-14,maxiters=50,verbose=true,exception=true,progress=false,
)

RB.solve!(
    prob,
    RB.ZhongQCCPNMono();
    dt,tspan,ftol=1e-14,maxiters=50,verbose=true,exception=true,progress=false,
)


plot_traj!(sc;showmesh=false,showground=false)
me = RB.mechanical_energy!(sc)
lines(me.E)
lines(me.V)
lines(me.T)