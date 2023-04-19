
function get_kinetic_energy(tg,sol)
    ke = VectorOfArray(Vector{Vector{Float64}}())
    rbs = tg.rigidbodies
    for (q,q̇) in zip(sol.qs,sol.q̇s)
        TR.reset_forces!(tg)
        TR.distribute_q_to_rbs!(tg,q,q̇)
        rbs_ke = [TR.kinetic_energy_coords(rb) for rb in rbs]
        push!(ke.u,rbs_ke)
    end
    ke
end

function deformation_to_elastic_potential_energy(ds,l0,u0,string)
    ls = l0 .+ ds
    epe = Vector{Float64}()
    for l in ls
        Δl = l - u0
        k = string.k
        string_epe = 1/2*k*Δl^2
        push!(epe,string_epe)
    end
    epe
end

function get_elastic_potential_energy(tg,sol)
    epe = VectorOfArray(Vector{Vector{Float64}}())
    for q in sol.qs
        TR.reset_forces!(tg)
        TR.distribute_q_to_rbs!(tg,q)
        TR.update_cables_apply_forces!(tg)
        epe_cables = TR.potential_energy.(tg.cables)
        push!(epe,epe_cables)
    end
    epe
end

function velocity_to_kinetic_energy(vs,ωs,rb)
    inertia = rb.prop.inertia
    mass = rb.prop.mass
    rb_ke = Vector{Float64}()
    for (v,ω) in zip(vs,ωs)
        trans_ke = 1/2*mass*v^2
        rot_ke = 1/2*inertia*ω^2
        ke = trans_ke + rot_ke
        push!(rb_ke,ke)
    end
    rb_ke
end

function position_to_gravity_potential_energy(rs,rb)
    pe = Vector{Float64}()
    g = TR.get_gravity(rb)
    for r in rs
        push!(pe,-r*g[2]*rb.prop.mass)
    end
    pe
end

function get_gravity_potential_energy(tg,sol)
    gpe = VectorOfArray(Vector{Vector{Float64}}())
    rbs = tg.rigidbodies[2:3]
    for (q,q̇) in zip(sol.qs,sol.q̇s)
        TR.reset_forces!(tg)
        TR.distribute_q_to_rbs!(tg,q,q̇)
        rbs_gpe = [TR.potential_energy(rb) for rb in rbs]
        push!(gpe.u,rbs_gpe)
    end
    gpe
end

function adams_to_energy(res,tr,l0,u0)
    @unpack tg = tr
    rbs = tg.rigidbodies
    ss = tg.cables
    rb2_gpe = position_to_gravity_potential_energy(res("H2"),rbs[2])
    rb3_gpe = position_to_gravity_potential_energy(res("H3"),rbs[3])
    gpe = rb2_gpe .+ rb3_gpe
    rb2_ke = velocity_to_kinetic_energy(res("V2"),deg2rad.(res("W2")),rbs[2])
    rb3_ke = velocity_to_kinetic_energy(res("V3"),deg2rad.(res("W3")),rbs[3])
    ke = rb2_ke .+ rb3_ke
    ss_epe = [deformation_to_elastic_potential_energy(res("deformation",i),l0[i],u0[i],ss[i]) for i = 1:4]
    epe = sum(ss_epe)
    energy =  gpe .+ ke .+ epe
    energy_error = abs.(energy.-energy[1])./energy[1]
    ke,epe,gpe,energy,energy_error
end

function get_energy(tg,sol)
    epe = get_elastic_potential_energy(tg,sol)
    ke = get_kinetic_energy(tg,sol)
    gpe = get_gravity_potential_energy(tg,sol)
    epes = sum.(epe.u)
    kes = sum.(ke.u)
    gpes = sum.(gpe.u)
    energy = epes .+ kes .+ gpes
    kes,epes,gpes,energy
end
