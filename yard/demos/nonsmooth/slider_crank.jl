figdir::String = ""
if Sys.iswindows()
    figdir::String = raw"C:\Users\luo22\OneDrive\Papers\FrictionalContact\CMAME"
elseif Sys.isapple()
    figdir::String = raw"/Users/jacob/Library/CloudStorage/OneDrive-SharedLibraries-onedrive/Papers/FrictionalContact/CMAME"
end
#-- point mass end

#--  Spinning top
includet("deps.jl")
import Rible as RB
include("../../vis.jl")
includet("../../vis.jl")
include("../../../examples/robots/slider_crank.jl")
includet("../../../examples/robots/slider_crank.jl")
sc = slider_crank(;coordsType=RB.QCF.QC)
plot_traj!(sc;showmesh=false,showground=false)
dt = 1e-4
tspan = (0.0,1.0)
prob = RB.SimProblem(sc,(bot)->RB.dynfuncs(bot;gravity=true))
RB.solve!(
    prob,
    RB.Zhong06();
    dt,tspan,ftol=1e-14,maxiters=50,verbose=true,exception=true,progress=false,
)
plot_traj!(sc;showmesh=false,showground=false)
me = RB.mechanical_energy!(sc)
lines(me.E)
lines(me.V)
lines(me.T)