abstract type AbstractFlexibleBody{N,T} <: AbstractBody{N,T}end
abstract type AbstractFlexibleBodyProperty{N,T} <: AbstractBodyProperty{N,T} end
abstract type AbstractFlexibleBodyState{N,T} <: AbstractBodyState{N,T} end

struct FlexibleBodyProperty{N,T} <: AbstractFlexibleBodyProperty{N,T}
    id::Int
    type::Symbol
    mass::T
    r̄g::SVector{N,T}
    r̄ps::Vector{SVector{N,T}}
    μs::Vector{T}
    es::Vector{T}
end

function FlexibleBodyProperty(
        id,type,mass::T,r̄g,r̄ps,
        μs=zeros(T,length(r̄ps)),
        es=zeros(T,length(r̄ps));
    ) where T
    FlexibleBodyProperty(id,type,mass,r̄g,r̄ps,μs,es)
end

struct FlexibleBodyCoordinatesCache{fT,MT,JT,VT,eT,ArrayT,N,M}
    pres_idx::Vector{Int}
    free_idx::Vector{Int}
    pres_idx_mask::SVector{N,Bool}
    free_idx_mask::SVector{N,Bool}
	Φ_mask::SVector{M,Bool}
    nΦ::Int
    funcs::fT
    e::eT
    ė::eT
    M::MT
	M⁻¹::MT
	∂Mė∂e::JT
	∂M⁻¹p∂e::JT
	Ṁė::VT
    ∂T∂eᵀ::VT
    So::ArrayT
    Sg::ArrayT
    Sps::Vector{ArrayT}
end

function get_Φ_mask(ancs::ANCF.ANC)
	nΦ = ANCF.get_nconstraints(ancs)
	@SVector ones(Bool,nΦ)
end

function get_idx_mask(ancs::ANCF.ANC,pres_idx)
    ne = ANCF.get_ncoords(ancs)
    pres_idx_mask = SVector{ne}(i in pres_idx for i = 1:ne)
    free_idx_mask = .!pres_idx_mask
    free_idx_mask, pres_idx_mask
end

function get_CoordinatesCache(
        prop::FlexibleBodyProperty{N,T},
        ancs::ANCF.ANC,
        e::AbstractVector,
        ė=zero(e);
        pres_idx=Int[],
        Φ_mask=get_Φ_mask(ancs)
    ) where {N,T}
    (;mass,r̄ps) = prop
    nr̄ps = length(r̄ps)
    nΦ = sum(Φ_mask)
    free_idx_mask, pres_idx_mask = get_idx_mask(ancs,pres_idx)
    pres_idx = findall(pres_idx_mask)
    free_idx = findall(free_idx_mask)
    cf = ANCF.CoordinateFunctions(ancs,free_idx_mask,Φ_mask)
    M = ANCF.build_M(ancs)
    M⁻¹ = inv(M)
    ∂Mq̇∂q = zero(M)
    ∂M⁻¹p∂q = zero(M)
    Ṁq̇ = @MVector zeros(T,size(M,2))
    ∂T∂qᵀ = @MVector zeros(T,size(M,2))
    (;S,x) = cf
    So = S(x(zeros(N)))
    Sg = ANCF.build_Sg(ancs,mass)
    Sps = [typeof(So)(S(x(r̄ps[i]))) for i in 1:nr̄ps]
    FlexibleBodyCoordinatesCache(
        pres_idx,
        free_idx,
        pres_idx_mask,
        free_idx_mask,
        Φ_mask,nΦ,
        cf,e,ė,
        M,M⁻¹,
        ∂Mq̇∂q,∂M⁻¹p∂q,Ṁq̇,
        ∂T∂qᵀ,So,Sg,Sps)
end

struct FlexibleBodyState{N,T,cacheType,contactType} <: AbstractFlexibleBodyState{N,T}
    rg::MVector{N,T}
	"Velocity of the centroid"
    ṙg::MVector{N,T}
	"Positions of anchor points"
    rps::Vector{MVector{N,T}}
	"Velocity of anchor Points"
    ṙps::Vector{MVector{N,T}}
	"Resultant force"
    f::MVector{N,T}
	"Action forces on anchor points"
    fps::Vector{MVector{N,T}}
	"Cache"
    cache::cacheType
    "Contacts"
    contacts::Vector{contactType}
end

function FlexibleBodyState(prop::FlexibleBodyProperty{N,T},
        ancs::ANCF.ANC,
        e::AbstractVector,
        ė=zero(e);
        pres_idx=Int[],
        Φ_mask=get_Φ_mask(ancs)) where {N,T}
    (;r̄g,r̄ps,μs,es) = prop
    nr̄ps = length(r̄ps)
    nc = ANCF.get_ncoords(ancs)
    cache = get_CoordinatesCache(
        prop,
        ancs,
        MVector{nc}(e),
        MVector{nc}(ė);
        pres_idx,Φ_mask
    )
    (;Sg,Sps) = cache
    rg = MVector{N}(Sg*e)
    ṙg = MVector{N}(Sg*ė)
    f = @MVector zeros(T,N)
    rps = [MVector{N}(Sps[i]*e) for i in 1:nr̄ps]
    ṙps = [MVector{N}(Sps[i]*ė) for i in 1:nr̄ps]
    fps = [@MVector zeros(T,N) for i in 1:nr̄ps]
    contacts = [Contact(i,μs[i],es[i]) for i in 1:nr̄ps]
    FlexibleBodyState(rg,ṙg,rps,ṙps,f,fps,cache,contacts)
end

struct FlexibleBody{N,T,cacheType,meshType} <: AbstractFlexibleBody{N,T}
	"Properties of a flexible body"
    prop::FlexibleBodyProperty{N,T}
	"State of a flexible body"
    state::FlexibleBodyState{N,T,cacheType}
    "Mesh"
    mesh::meshType
end

function FlexibleBody(prop,state)
    FlexibleBody(prop,state,nothing)
end


function body2coordinates(fb::FlexibleBody)
    (;e,ė) = fb.state.cache
    e,ė
end

function update_rigid!(
        state::FlexibleBodyState,
        cache::FlexibleBodyCoordinatesCache{<:ANCF.CoordinateFunctions},
        prop::AbstractBodyProperty,q,q̇)
    (;cache,rg,ṙg) = state
    (;e, ė, Sg) = cache
    e .= q
    ė .= q̇
    mul!(rg, Sg, q)
    mul!(ṙg, Sg, q̇)
end

function stretch_rigid!(
    cache::FlexibleBodyCoordinatesCache,
	prop::AbstractBodyProperty,c)
end

function update_transformations!(
    cache::FlexibleBodyCoordinatesCache,
    state::FlexibleBodyState,
    prop::AbstractBodyProperty,q)
end

function move_rigid!(
        state::FlexibleBodyState,
        cache::FlexibleBodyCoordinatesCache{<:ANCF.CoordinateFunctions},
        prop::AbstractBodyProperty,e,ė)
    (;rps,ṙps) = state
    (;Sps) = cache
    for (i,(rp,ṙp)) in enumerate(zip(rps,ṙps))
        mul!(rp, Sps[i], e)
        mul!(ṙp, Sps[i], ė)
    end
end

function generalize_force!(F,state::AbstractFlexibleBodyState)
	(;cache,f,fps) = state
	(;Sps,Sg,e,∂T∂eᵀ) = cache
	for (pid,fp) in enumerate(fps)
		# F .+= transpose(Sps[pid])*fp
		mul!(F,transpose(Sps[pid]),fp,1,1)
	end
	# F .+= transpose(Sg)*f
	mul!(F,transpose(Sg),f,1,1)
	F .+= ∂T∂eᵀ
    Q = ANCF.make_Q(cache.funcs.ancs)
	F .-= Q(e)
	F
end

"""
Return the potential energy for gravity of a flexible body
$(TYPEDSIGNATURES)
"""
function potential_energy_gravity(fb::AbstractFlexibleBody)
	(;mass) = fb.prop
    (;rg) = fb.state
    gravity_acceleration = get_gravity(fb)
    -transpose(rg)*gravity_acceleration*mass
end

"""
Return the potential energy for strain of a flexible body
$(TYPEDSIGNATURES)
"""
function potential_energy_strain(fb::AbstractFlexibleBody)
    (;cache) = fb.state
    (;funcs,e) = cache
    (;ancs) = funcs
    ANCF.make_V(ancs)(e)
end

function subdivide(fb::FlexibleBody,nx=1,ny=1,nz=1)
    (;prop,state) = fb
    (;cache) = state
    (;pres_idx,funcs,e,ė) = cache
    (;ancs) = funcs
    (;radius,E,L,ρ) = ancs
    xrange = range(0,L,nx+1) |> collect
    _S = ANCF.make_S(ancs)
    _Sₓ = ANCF.make_Sₓ(ancs)
    rgrid = _S.(xrange).*Ref(e)
    rₓgrid = _Sₓ.(xrange).*Ref(e)
    ṙgrid = _S.(xrange).*Ref(ė)
    ṙₓgrid = _Sₓ.(xrange).*Ref(ė)
    esegs = [
        [rgrid[1];rₓgrid[1];rgrid[2];rₓgrid[2]]
    ]
    append!(esegs,[
        [rgrid[i];rₓgrid[i];rgrid[i+1];rₓgrid[i+1]]
        for i = 2:nx
    ])
    ėsegs = [
        [ṙgrid[1];ṙₓgrid[1];ṙgrid[2];ṙₓgrid[2]]
    ]
    append!(ėsegs,[
        [ṙgrid[i];ṙₓgrid[i];ṙgrid[i+1];ṙₓgrid[i+1]]
        for i = 2:nx
    ])
    T = typeof(L)
    # @show pres_idx
    subfbs = [
        begin            
            subL = xrange[i+1] - xrange[i]
            subancs = ANCF.ANC3DRURU(ρ;E,L=subL,radius)
            mass = ANCF.build_mass(subancs,(x)->ρ(xrange[i+1]+x))
            r̄g = SVector(subL/2,T(0),T(0))
            r̄p1 = SVector(T(0),T(0),T(0))
            r̄p2 = SVector(subL,T(0),T(0))
            r̄ps = [
                r̄p1,r̄p2
            ]
            prop = FlexibleBodyProperty(
                i,
                :cable,
                mass,
                r̄g,
                # length(r̄ps),
                r̄ps
            )
            if i == 1
                sub_pres_idx = filter((j)->j in 1:6, pres_idx)
            elseif i == nx
                sub_pres_idx = filter((j)->j in 7:12, pres_idx)
            else
                sub_pres_idx = filter((j)->false, pres_idx)
            end
            # @show sub_pres_idx
            # cache = get_CoordinatesCache(prop,ancs,𝐞)
            state = FlexibleBodyState(prop,subancs,esegs[i],ėsegs[i];pres_idx=sub_pres_idx)
            fb = FlexibleBody(prop,state)
        end
        for i = 1:nx
    ]
    sm = reduce(
        vcat,[
            begin
                ssm = zeros(Int,6,nx)
                ssm[1:6,i] = 7:12
                ssm[1:6,i+1] = 1:6
                ssm
            end
            for i = 1:nx-1
        ]
    )
    subfbs,sm
end

function clear_forces!(fb::AbstractFlexibleBody)
    (;state) = fb
    state.f .= 0
    foreach(state.fps) do fp
        fp .= 0
    end
end
