function dynfuncs1(bot)
    (;tg) = bot
    function F!(F,q,q̇,s,t)
        TR.clear_forces!(tg)
        TR.update_rigids!(tg,q,q̇)
        TR.distribute_s̄!(tg,s)
        TR.update_tensiles!(tg)
        # TR.apply_gravity!(tg)
        TR.generate_forces!(tg)
        TR.get_force!(F,tg)
    end
    Jac_F! = true
    @eponymtuple(F!, Jac_F!)
end

function actuate1!(tg, t; dt=1e-2)
    function inner1(t;dt=dt)
        a = 3000.0
        if 0<t<5
            return -a*dt
        # elseif 1<=t<5
        #     return -a*dt
        # elseif 5<=t<6
        #     return a*t*dt - 6a*dt
        else
            return 0
        end
    end
    function inner2(t; dt = dt)
        a = 2000.0
        if 4<t<8
            return -a*dt
        else
            return 0
        end
    end
    # tg.tensiles.clustercables[1].segs[1].state.restlen -= inner1(t)
    tg.tensiles.clustercables[2].segs[1].state.restlen += inner1(t)
    tg.tensiles.clustercables[3].segs[1].state.restlen += inner1(t)
    # tg.tensiles.clustercables[5].segs[1].state.restlen += inner1(t)
    # tg.tensiles.clustercables[6].segs[1].state.restlen += inner1(t)
end


function actuate2!(tg, t; dt=1e-2)
    function inner1(t;dt=dt)
        a = 2550.0
        if 0<t<8
            return -a*dt
        # elseif 1<=t<5
        #     return -a*dt
        # elseif 5<=t<6
        #     return a*t*dt - 6a*dt
        else
            return 0
        end
    end
    function inner2(t; dt = dt)
        a = 3000.0
        if 3<t<7
            return -a*dt
        else
            return 0
        end
    end
    tg.tensiles.clustercables[1].segs[1].state.restlen += inner1(t)
    tg.tensiles.clustercables[2].segs[1].state.restlen -= inner1(t)
    tg.tensiles.clustercables[3].segs[1].state.restlen -= inner1(t)
end

function actuate2_1!(tg, t; dt=1e-2)
    function inner1(t;dt=dt)
        a = 1000
        if 0<t<10
            return -a*dt
        # elseif 1<=t<5
        #     return -a*dt
        # elseif 5<=t<6
        #     return a*t*dt - 6a*dt
        else
            return -a*dt
        end
    end
    function inner2(t; dt = dt)
        a = 1000.0
        if 3<t<7
            return -a*dt
        else
            return 0
        end
    end
    # tg.tensiles.clustercables[1].segs[1].state.restlen += 1inner1(t)
    tg.tensiles.clustercables[4].segs[1].state.restlen += 1inner1(t)
    tg.tensiles.clustercables[7].segs[2].state.restlen += .05inner1(t)
    tg.tensiles.clustercables[7].segs[3].state.restlen += .05inner1(t)
end


function actuate3!(tg, t; dt=1e-2)
    function inner1(t;dt=dt)
        a = 1000.0
        if 0<t<4
            return -a*dt
        else
            return 0
        end
    end
    function inner2(t; dt = dt)
        a = 1000.0
        if 5<t<9
            return -a*dt
        else
            return 0
        end
    end
    tg.tensiles.clustercables[7].segs[1].state.restlen += 2inner1(t)
    # tg.tensiles.clustercables[2].segs[2].state.restlen -= 0.2inner1(t)
    # tg.tensiles.clustercables[3].segs[2].state.restlen -= 0.2inner1(t)
    tg.tensiles.clustercables[2].segs[10].state.restlen += 3.0inner2(t)
    # tg.tensiles.clustercables[4].segs[1].state.restlen -= 0.3inner2(t)
    tg.tensiles.clustercables[5].segs[1].state.restlen += 1.8inner2(t)
    tg.tensiles.clustercables[6].segs[1].state.restlen += 2.0inner2(t)
end

function actuate_for_2!(tg, t; dt=1e-2)
    function inner1(t;dt=dt)
        a = 1000.0
        if 0<t<4
            return -a*dt
        else
            return 0
        end
    end
    tg.tensiles.clustercables[2].segs[1].state.restlen -= .5inner1(t)
    tg.tensiles.clustercables[3].segs[1].state.restlen += .5inner1(t)
    tg.tensiles.clustercables[1].segs[1].state.restlen += .5inner1(t)

end
function actuate_for_3!(tg, t; dt=1e-2)
    function inner1(t;dt=dt)
        a = 1000.0
        if 0<t<4
            return -a*dt
        else
            return 0
        end
    end
    tg.tensiles.clustercables[3].segs[1].state.restlen -= .5inner1(t)
    tg.tensiles.clustercables[2].segs[1].state.restlen += .5inner1(t)
    tg.tensiles.clustercables[1].segs[1].state.restlen += .5inner1(t)
end

function actuate_3!(tg, t; dt=1e-2)
    function inner1(t;dt=dt)
        a = 1000.0
        if 0<t<4
            return -a*dt
        else
            return 0
        end
    end
    tg.tensiles.clustercables[3].segs[1].state.restlen += 2.5inner1(t)
end
function actuate_2!(tg, t; dt=1e-2)
    function inner1(t;dt=dt)
        a = 1000.0
        if 0<t<4
            return -a*dt
        else
            return 0
        end
    end
    tg.tensiles.clustercables[2].segs[1].state.restlen += 2.5inner1(t)
end


function actuate11!(tg, t; dt=1e-2)
    function inner1(t;dt=dt)
        a = 3000.0
        if 0<t<5
            return -a*dt
        # elseif 1<=t<5
        #     return -a*dt
        # elseif 5<=t<6
        #     return a*t*dt - 6a*dt
        else
            return 0
        end
    end
    # tg.tensiles.clustercables[1].segs[1].state.restlen -= inner1(t)
    tg.tensiles.clustercables[1].segs[1].state.restlen += inner1(t)
end


function actuate12!(tg, t; dt=1e-2)
    function inner1(t;dt=dt)
        a = 1000.0
        if 0<t<3.2
            return -a*dt
        # elseif 1<=t<5
        #     return -a*dt
        # elseif 5<=t<6
        #     return a*t*dt - 6a*dt
        else
            return 0
        end
    end
    # tg.tensiles.clustercables[1].segs[1].state.restlen -= inner1(t)
    # tg.tensiles.clustercables[5].segs[1].state.restlen += 0.5inner1(t)
    # tg.tensiles.clustercables[6].segs[1].state.restlen += 0.5inner1(t)
    tg.tensiles.clustercables[7].segs[4].state.restlen += 1.2inner1(t)
end

function actuate121!(tg, t; dt=1e-2)
    function inner1(t;dt=dt)
        a = 3000.0
        if 0<t<5
            return -a*dt
        # elseif 1<=t<5
        #     return -a*dt
        # elseif 5<=t<6
        #     return a*t*dt - 6a*dt
        else
            return 0
        end
    end
    function inner2(t; dt = dt)
        a = 2000.0
        if 4<t<8
            return -a*dt
        else
            return 0
        end
    end
    # tg.tensiles.clustercables[1].segs[1].state.restlen -= inner1(t)
    tg.tensiles.clustercables[5].segs[1].state.restlen += .5inner1(t)
    tg.tensiles.clustercables[6].segs[1].state.restlen += .5inner1(t)
    tg.tensiles.clustercables[8].segs[1].state.restlen += .3inner1(t)
    tg.tensiles.clustercables[9].segs[1].state.restlen += .3inner1(t)
    # tg.tensiles.clustercables[5].segs[1].state.restlen += inner1(t)
    # tg.tensiles.clustercables[6].segs[1].state.restlen += inner1(t)
end