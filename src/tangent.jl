
function build_tangent(tg)
    q, q̇ = get_q(tg)
    J = [build_Ji(tg,i) for i = 1:tg.ncables]
    U = [transpose(Ji)*Ji for Ji in J]
    l = [sqrt(transpose(q)*Ui*q) for Ui in U]
    # display(l)
    # display([s.state.length for s in tg.cables])
    # @show l.-[s.state.length for s in tg.cables]
    l̇ = [(transpose(q)*Ui*q̇)/li for (Ui,li) in zip(U,l)]
    ∂l∂q = [(transpose(q)*Ui)./li for (Ui,li) in zip(U,l)]
    ∂l̇∂q̇ = ∂l∂q
    ∂l̇∂q = [(li*transpose(q̇)-l̇i*transpose(q))/li^2*Ui for (Ui,li,l̇i) in zip(U,l,l̇)]
    k = [s.k for s in tg.cables]
    c = [s.c for s in tg.cables]
    ∂f∂q = [ki*∂li∂q+ci*∂l̇i∂q for (ki,ci,∂li∂q,∂l̇i∂q) in zip(k,c,∂l∂q,∂l̇∂q)]
    ∂f∂q̇ = [ci*∂l̇i∂q̇ for (ci,∂l̇i∂q̇) in zip(c,∂l̇∂q̇)]
    l̂ = [Ji*q/li for (Ji,li) in zip(J,l)]
    ∂l̂∂q = [(Ji-l̂i*∂li∂q)/li for (li,l̂i,Ji,∂li∂q) in zip(l,l̂,J,∂l∂q)]
    # @show l̂.-[s.state.direction for s in tg.cables]
    f = [s.state.tension for s in tg.cables]
    ∂𝐟∂q = [l̂i*∂fi∂q+fi*∂l̂i∂q for (l̂i,fi,∂fi∂q,∂l̂i∂q) in zip(l̂,f,∂f∂q,∂l̂∂q)]
    ∂𝐟∂q̇ = [l̂i*∂fi∂q̇ for (l̂i,∂fi∂q̇) in zip(l̂,∂f∂q̇)]
    @unpack ndim,ncoords,ncables = tg
    ∂Γ∂q = zeros(eltype(q),ndim*ncables,ncoords)
    ∂Γ∂q̇ = zeros(eltype(q),ndim*ncables,ncoords)
    for i in 1:ncables
        ∂Γ∂q[(i-1)*ndim+1:i*ndim,:] = ∂𝐟∂q[i]
        ∂Γ∂q̇[(i-1)*ndim+1:i*ndim,:] = ∂𝐟∂q̇[i]
    end
    ∂Γ∂q,∂Γ∂q̇
end

function build_Jac_Γ(tg)
    ns = tg.ncables
    @unpack ncoords,ndim = tg
    J = [build_Ji(tg,i) for i = 1:ns]
    U = [transpose(Ji)*Ji for Ji in J]
    k = [s.k for s in tg.cables]
    c = [s.c for s in tg.cables]
    function inner_Jac_Γ(q,q̇)
        reset_forces!(tg)
        distribute_q_to_rbs!(tg,q,q̇)
        update_cables_apply_forces!(tg)
        f = [s.state.tension for s in tg.cables]
        l = [s.state.length for s in tg.cables]
        u = [s.state.restlen for s in tg.cables]
        qᵀU = [transpose(q)*U[i] for i = 1:ns]
        # l = [sqrt(qᵀU[i]*q) for i = 1:ns]
        l̇ = [(qᵀU[i]*q̇)/l[i] for i = 1:ns]
        l̂ = [J[i]*q/l[i] for i = 1:ns]
        ∂Γ∂q = zeros(eltype(q),ndim*ns,ncoords)
        ∂Γ∂q̇ = zeros(eltype(q),ndim*ns,ncoords)
        for i in 1:ns
            l̂qᵀUi = l̂[i]*qᵀU[i]
            # ∂Γ∂q[(i-1)*ndim+1:i*ndim,:] .= l̂[i]*(k[i]./l[i].*qᵀU[i] .+ c[i]./l[i].*transpose(q̇).-c[i]*l̇[i]./l[i]^2 .*qᵀU[i])
            # ∂Γ∂q[(i-1)*ndim+1:i*ndim,:] .+= f[i].*(J[i]./l[i].-l̂[i]*qᵀU[i]./l[i]^2)
            ∂Γ∂q[(i-1)*ndim+1:i*ndim,:] .= (k[i]*u[i]-c[i]*l̇[i])/l[i]^2 .*l̂qᵀUi .+ c[i]/l[i].*(l̂[i]*transpose(q̇)*U[i]) .+ f[i]/l[i] .*J[i]
            ∂Γ∂q̇[(i-1)*ndim+1:i*ndim,:] .= c[i]/l[i].*l̂qᵀUi
        end
        ∂Γ∂q,∂Γ∂q̇
    end
end

function make_testtangent(tgstruct)
    @unpack ncoords,ncables = tgstruct
    function build_𝐟(x)
        q = x[1:ncoords]
        q̇ = x[ncoords+1:2ncoords]
        reset_forces!(tgstruct)
        distribute_q_to_rbs!(tgstruct,q,q̇)
        update_cables_apply_forces!(tgstruct)
        vcat([tgstruct.cables[i].state.direction*tgstruct.cables[i].state.tension
            for i = 1:ncables]...)
    end
end
