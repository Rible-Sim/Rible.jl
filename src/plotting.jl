function get_strings(tgstruct)
    string2ap = tgstruct.connectivity.string2ap
    rbs = tgstruct.rigidbodies
    [Point(rbs[s[1].rbid].state.rps[s[1].apid]) =>
     Point(rbs[s[2].rbid].state.rps[s[2].apid])
     for s in string2ap]
end

function get_cables(tgstruct)
    string2ap = tgstruct.connectivity.string2ap
    rbs = tgstruct.rigidbodies
    cables = [Point(rbs[s[1].rbid].state.rps[s[1].apid]) =>
     Point(rbs[s[2].rbid].state.rps[s[2].apid])
     for s in string2ap]
end

function get_clustercables(tgstruct)
    clusterstring2ap = Vector{Tuple{TensegrityRobots.ID, TensegrityRobots.ID}}()
    temp = tgstruct.connectivity.clusterstring2ap
    for i in 1:length(temp)
        append!(clusterstring2ap, temp[i])
    end
    rbs = tgstruct.rigidbodies
    clustercables = [Point(rbs[s[1].rbid].state.rps[s[1].apid]) =>
     Point(rbs[s[2].rbid].state.rps[s[2].apid])
     for s in clusterstring2ap]
end