get_strings(tg) = get_cables(tg)

function get_cables(tg)
    string2ap = tg.connectivity.string2ap
    rbs = tg.rigidbodies
    cables = [Point(rbs[s[1].rbid].state.rps[s[1].apid]) =>
     Point(rbs[s[2].rbid].state.rps[s[2].apid])
     for s in string2ap]
end

function get_clustercables(tg)
    clusterstring2ap = tg.connectivity.clusterstring2ap
    rbs = tg.rigidbodies
    clustercables = [
    [
    Point(rbs[s[1].rbid].state.rps[s[1].apid]) => Point(rbs[s[2].rbid].state.rps[s[2].apid]) 
    for s in cs
    ] for cs in clusterstring2ap
    ]
end
