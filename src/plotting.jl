function get_strings(tg::TensegrityStructure)
    string2ap = tg.connectivity.string2ap
	T = get_numbertype(tg)
	ndim = get_ndim(tg)
    if ndim == 2
		ret = Vector{Pair{Point2{T},Point2{T}}}(undef,tg.nstrings)
	elseif ndim == 3
		ret = Vector{Pair{Point3{T},Point3{T}}}(undef,tg.nstrings)
	else
		error("Somethings no right")
	end
    map!(ret,string2ap) do scnt
        Point(scnt.end1.rbsig.state.rps[scnt.end1.pid]) =>
        Point(scnt.end2.rbsig.state.rps[scnt.end2.pid])
    end
    ret
end


function make_apcnt(rigidbodies)
	ret = Vector{Vector{Int}}()
	is = Ref(0)
	foreach(rigidbodies) do rb
		push!(ret,collect(is[]+1:is[]+rb.prop.naps))
		is[] += rb.prop.naps
	end
	ret
end
