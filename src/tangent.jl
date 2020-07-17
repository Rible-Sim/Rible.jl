struct TensionTangent{T}
    J::Vector{Array{T,2}}
    U::Vector{Array{T,2}}
    l::Vector{T}
    l̇::Vector{T}
    ∂l∂q::Vector{Array{T,2}}
    ∂l̇∂q::Vector{Array{T,2}}
    ∂l̇∂q̇::Vector{Array{T,2}}
    ∂f∂q::Vector{Array{T,2}}
    ∂f∂q̇::Vector{Array{T,2}}
    ∂l̂∂q::Vector{Array{T,2}}
end

function build_tangent(tgstruct,q,q̇)
    reset_forces!(tgstruct)
    distribute_q_to_rbs!(tgstruct,q,q̇)
    update_strings_apply_forces!(tgstruct)
    J = [build_Ji(tgstruct,i) for i = 1:tgstruct.nstrings]
    U = [transpose(Ji)*Ji for Ji in J]
    l = [sqrt(transpose(q)*Ui*q) for Ui in U]
    # display(l)
    # display([s.state.length for s in tgstruct.strings])
    # @show l.-[s.state.length for s in tgstruct.strings]
    l̇ = [(transpose(q)*Ui*q̇)/li for (Ui,li) in zip(U,l)]
    ∂l∂q = [(transpose(q)*Ui)./li for (Ui,li) in zip(U,l)]
    ∂l̇∂q̇ = ∂l∂q
    ∂l̇∂q = [(li*transpose(q̇)-l̇i*transpose(q))/li^2*Ui for (Ui,li,l̇i) in zip(U,l,l̇)]
    k = [s.k for s in tgstruct.strings]
    c = [s.c for s in tgstruct.strings]
    ∂f∂q = [ki*∂li∂q+ci*∂l̇i∂q for (ki,ci,∂li∂q,∂l̇i∂q) in zip(k,c,∂l∂q,∂l̇∂q)]
    ∂f∂q̇ = [ci*∂l̇i∂q̇ for (ci,∂l̇i∂q̇) in zip(c,∂l̇∂q̇)]
    l̂ = [Ji*q/li for (Ji,li) in zip(J,l)]
    ∂l̂∂q = [(Ji-l̂i*∂li∂q)/li for (li,l̂i,Ji,∂li∂q) in zip(l,l̂,J,∂l∂q)]
    # @show l̂.-[s.state.direction for s in tgstruct.strings]
    f = [s.state.tension for s in tgstruct.strings]
    ∂𝐟∂q = [l̂i*∂fi∂q+fi*∂l̂i∂q for (l̂i,fi,∂fi∂q,∂l̂i∂q) in zip(l̂,f,∂f∂q,∂l̂∂q)]
    ∂𝐟∂q̇ = [l̂i*∂fi∂q̇ for (l̂i,∂fi∂q̇) in zip(l̂,∂f∂q̇)]
    @unpack ndim,ncoords,nstrings = tgstruct
    ∂L∂q = zeros(eltype(q),ndim*nstrings,ncoords)
    ∂L∂q̇ = zeros(eltype(q),ndim*nstrings,ncoords)
    for i in 1:nstrings
        ∂L∂q[(i-1)*ndim+1:i*ndim,:] = ∂𝐟∂q[i]
        ∂L∂q̇[(i-1)*ndim+1:i*ndim,:] = ∂𝐟∂q̇[i]
    end
    ∂L∂q,∂L∂q̇
end

function make_testtangent(tgstruct)
    @unpack ncoords,nstrings = tgstruct
    function build_𝐟(x)
        q = x[1:ncoords]
        q̇ = x[ncoords+1:2ncoords]
        reset_forces!(tgstruct)
        distribute_q_to_rbs!(tgstruct,q,q̇)
        update_strings_apply_forces!(tgstruct)
        vcat([tgstruct.strings[i].state.direction*tgstruct.strings[i].state.tension
            for i = 1:nstrings]...)
    end
end
