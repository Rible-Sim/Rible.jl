
function make_spine(n,θ=0.0)
    a = 0.1
    b = 0.1*√2
    α = 3/4*π
    d = 0.15

    function vert(i,r,θ,a,b,α)
        if i == 1
            contactable = false
            visible = true
            ci = collect(1:6)
			cstr_idx = Int[]
        else
            contactable = true
            visible = true
            ci = Int[]
			cstr_idx = collect(1:3)
	        # ap1 = b*[cos( α),sin( α)]
	        # ap2 = b*[cos( α),sin( α)]
	        # ap3 = b*[cos(-α),sin(-α)]
	        # ap4 = b*[cos(-α),sin(-α)]
        end		
		ap1 = [a,0.0]
		ap2 = b*[cos( α),sin( α)]
		ap3 = b*[cos(-α),sin(-α)]
		# ap4 = [a,0.0]
        mass_locus = zeros(2)
        aps = [ap1,ap2,ap3]
        m = 0.495 #kg
        inertia = SMatrix{2,2}(Diagonal(ones(2)))
        prop = RB.RigidBodyProperty(
			i,
			contactable,m,inertia,mass_locus,aps;
			visible
		)

        ri = SVector{2}(r)
        ro = SVector{2}(r)
        ω = 0.0
        ṙo = @SVector zeros(2)
        nmcs,_ = RB.NCF.NC1P2V(ri,ro,θ)
        state = RB.RigidBodyState(prop,ro,θ,ṙo,ω, ci, cstr_idx)
        RB.RigidBody(prop,state)
    end

    rs = [[i*d,0.0] for i = 0:n-1]
    θs = [i*θ for i = 0:n-1]
    rbs = [vert(i,rs[i],θs[i],a,b,α) for i = 1:n]

	rigdibodies = TypeSortedCollection(rbs)
    numbered = RB.number(rigdibodies)
    indexed = RB.index(rigdibodies)

    ncables = 4*(n-1)
    k = 840.0 #N/m
    c = 100.0
    original_restlens = repeat([0.15, 0.11, 0.11, 0.15]./2,n-1)
    ss = [RB.DistanceSpringDamper2D(original_restlens[i],k,c) for i = 1:ncables]
    apparatuses = (cables = ss,)
    hub = nothing
	

	matrix_cnt = zeros(Int,ncables,n)
	for i = 1:n-1
		matrix_cnt[4(i-1)+1,[i,i+1]] = [2,-2]
		matrix_cnt[4(i-1)+2,[i,i+1]] = [3,-3]
		matrix_cnt[4(i-1)+3,[i,i+1]] = [1,-2]
		matrix_cnt[4(i-1)+4,[i,i+1]] = [1,-3]
	end
	matrix_cnt |> display
	connected = RB.connect(rigdibodies, matrix_cnt)
	tensioned = @eponymtuple(connected,)

    cnt = RB.Connectivity(numbered,indexed,tensioned)
	# rbs
    st = RB.Structure(rigdibodies,apparatuses,cnt)
    bot = RB.Robot(st,hub)
end
