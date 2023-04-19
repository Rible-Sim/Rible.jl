using Setfield
# Setfield.jl is used to provide the parameter axis @lens
import BifurcationKit as BK

N(x; a = 0.5, b = 0.01) = 1 + (x + a*x^2)/(1 + b*x^2)
function F_chan(x, p)
	(; α, β) = p
	f = similar(x)
	n = length(x)
	f[1] = x[1] - β
	f[n] = x[n] - β
	for i=2:n-1
		f[i] = (x[i-1] - 2 * x[i] + x[i+1]) * (n-1)^2 + α * N(x[i], b = β)
	end
	return f
end
n = 101
sol0 = [(i-1)*(n-i)/n^2+0.1 for i=1:n]

# set of parameters
par = (α = 3.3, β = 0.01)

optnewton = BK.NewtonPar(tol = 1e-11, verbose = true)

prob = BK.BifurcationProblem(F_chan, sol0, par, (@lens _.α),)
prob.VF
sol = @time BK.newton( prob, optnewton)
optcont = BK.ContinuationPar(dsmin = 0.01, dsmax = 0.2, ds= 0.1, pMin = 0., pMax = 4.2,
newtonOptions = BK.NewtonPar(maxIter = 10, tol = 1e-9))

br = BK.continuation(prob, BK.PALC(), optcont; plot = false)
# index of the Fold bifurcation point in br.specialpoint
indfold = 2

outfold = BK.newton(br, indfold)

BK.converged(outfold)
printstyled(color=:red, "--> We found a Fold Point at α = ", outfold.u.p, ", β = 0.01, from ", br.specialpoint[indfold].param,"\n")

outfoldco = BK.continuation(br, indfold,
	# second parameter axis to use for codim 2 curve
	(@lens _.β),
	# we disable the computation of eigenvalues, it makes little sense here
	BK.ContinuationPar(optcont, detectBifurcation = 0))


scene = plot(outfoldco, plotfold = true, legend = :bottomright)
