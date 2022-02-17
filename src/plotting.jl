function get_strings(tg::TensegrityStructure{T,N}) where {T,N}
    string2ap = tg.connectivity.string2ap
    ret = Vector{Pair{Point2{Float64},Point2{Float64}}}(undef,tg.nstrings)
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
