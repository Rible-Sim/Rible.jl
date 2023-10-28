
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

function get_elastic_potential_energy(st,sol)
    epe = VectorOfArray(Vector{Vector{Float64}}())
    for q in sol.qs
        RB.reset_forces!(st)
        RB.distribute_q_to_rbs!(st,q)
        RB.update_strings_apply_forces!(st)
        epe_strings = RB.potential_energy.(st.strings)
        push!(epe,epe_strings)
    end
    epe
end

function velocity_to_kinetic_energy(vs,ω1s,ω2s,ω3s,rb)
    inertia = body.prop.inertia
    mass = body.prop.mass
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

function get_kinetic_energy(st,sol)
    ke = VectorOfArray(Vector{Vector{Float64}}())
    rbs = st.rigidbodies
    for (q,q̇) in zip(sol.qs,sol.q̇s)
        RB.reset_forces!(st)
        RB.distribute_q_to_rbs!(st,q,q̇)
        rbs_ke = [RB.kinetic_energy_coords(body) for rb in rbs]
        push!(ke.u,rbs_ke)
    end
    ke
end

function cm_to_gravity_potential_energy(cms,rb)
    g = RB.get_gravity(body)
    gpes = -g[end]*cms*body.prop.mass
end

function adams_to_energy(res,tr,l0,s_index,energy1;gravity=false)
    @unpack st = tr
    rbs = st.rigidbodies
    ss = st.strings
    u0 = RB.get_original_restlen(tr)
    ke = velocity_to_kinetic_energy(res("PV"),deg2rad.(res("PWX")),deg2rad.(res("PWY")),deg2rad.(res("PWZ")),st.rigidbodies[2])
    ss_epe = [deformation_to_elastic_potential_energy(res("deformation",s_index[i]),l0[i],u0[i],ss[i]) for i = 1:6]
    epe = sum(ss_epe)
    energy =  ke .+ epe
    if gravity
        gpe = cm_to_gravity_potential_energy(res("P2_CM"),st.rigidbodies[2])
        energy .+= gpe
    end
    energy1 = energy[1]
    energy_err = abs.((energy.-energy1)./energy1)
    ke,epe,energy,energy_err
end
