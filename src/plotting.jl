function get_strings(tgstruct)
    string2ap = tgstruct.connectivity.string2ap
    rbs = tgstruct.rigidbodies
    [Point(rbs[s[1].rbid].state.p[s[1].apid]) =>
     Point(rbs[s[2].rbid].state.p[s[2].apid])
     for s in string2ap]
end