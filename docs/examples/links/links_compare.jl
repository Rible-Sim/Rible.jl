
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
        TR.update_strings_apply_forces!(tg)
        epe_strings = TR.potential_energy.(tg.strings)
        push!(epe,epe_strings)
    end
    epe
end

function velocity_to_kinetic_energy(vs,ω1s,ω2s,ω3s,rb)
    inertia = rb.prop.inertia
    mass = rb.prop.mass
    rb_ke = Vector{Float64}()
    for (v,ω1,ω2,ω3) in zip(vs,ω1s,ω2s,ω3s)
        trans_ke = 1/2*mass*v^2
        ω = [ω1,ω2,ω3]
        rot_ke = 1/2*transpose(ω)*inertia*ω
        ke = trans_ke + rot_ke
        push!(rb_ke,ke)
    end
    rb_ke
end

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

function cm_to_gravity_potential_energy(cms,rb)
    g = TR.get_gravity(rb)
    gpes = -g[end]*cms*rb.prop.mass
end

function adams_to_energy(res,tr,l0,s_index,energy1;gravity=false)
    @unpack tg = tr
    rbs = tg.rigidbodies
    ss = tg.strings
    u0 = TR.get_original_restlen(tr)
    ke = velocity_to_kinetic_energy(res("PV"),deg2rad.(res("PWX")),deg2rad.(res("PWY")),deg2rad.(res("PWZ")),tg.rigidbodies[2])
    ss_epe = [deformation_to_elastic_potential_energy(res("deformation",s_index[i]),l0[i],u0[i],ss[i]) for i = 1:6]
    epe = sum(ss_epe)
    energy =  ke .+ epe
    if gravity
        gpe = cm_to_gravity_potential_energy(res("P2_CM"),tg.rigidbodies[2])
        energy .+= gpe
    end
    energy1 = energy[1]
    energy_err = abs.((energy.-energy1)./energy1)
    ke,epe,energy,energy_err
end
