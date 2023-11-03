abstract type AbstractFlexibleBody{N,T} <: AbstractBody{N,T}end
abstract type AbstractFlexibleBodyProperty{N,T} <: AbstractBodyProperty{N,T} end
abstract type AbstractFlexibleBodyState{N,T} <: AbstractBodyState{N,T} end

struct FlexibleBodyProperty{N,T} <: AbstractFlexibleBodyProperty{N,T}
    id::Int
    type::Symbol
    mass::T
    mass_locus::Locus{N,T}
    loci::Vector{Locus{N,T}}
end

function FlexibleBodyProperty(
        id,type,mass::T,
        mass_center_position::SVector{N,T},positions,
        friction_coefficients=zeros(T,length(positions)),
        restitution_coefficients=zeros(T,length(positions));
    ) where {N,T}
    VecType = SVector{N}
    mass_locus = Locus(
        VecType(mass_center_position)
    )
    loci = [
        Locus(
            VecType(position),
            friction_coefficient,
            restitution_coefficient
        )
        for (position,friction_coefficient,restitution_coefficient) in zip(
            positions,
            friction_coefficients,
            restitution_coefficients
        )
    ]
    FlexibleBodyProperty(
        id,type,
        mass,
        mass_locus,
        loci
    )
end

struct FlexibleBodyCoordinatesCache{fT,MT,JT,VT,eT,ArrayT,N,M}
    pres_idx::Vector{Int}
    free_idx::Vector{Int}
    pres_idx_mask::SVector{N,Bool}
    free_idx_mask::SVector{N,Bool}
	cstr_idx::SVector{M,Bool}
    num_of_cstr::Int
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

function get_cstr_idx(ancs::ANCF.ANC)
	num_of_cstr = ANCF.get_num_of_cstr(ancs)
	@SVector ones(Bool,num_of_cstr)
end

function get_idx_mask(ancs::ANCF.ANC,pres_idx)
    ne = ANCF.get_num_of_coords(ancs)
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
        cstr_idx=get_cstr_idx(ancs)
    ) where {N,T}
    (;mass,loci) = prop
    num_of_loci = length(loci)
    num_of_cstr = sum(cstr_idx)
    free_idx_mask, pres_idx_mask = get_idx_mask(ancs,pres_idx)
    pres_idx = findall(pres_idx_mask)
    free_idx = findall(free_idx_mask)
    cf = ANCF.CoordinateFunctions(ancs,free_idx_mask,cstr_idx)
    M = ANCF.build_M(ancs)
    M⁻¹ = inv(M)
    ∂Mq̇∂q = zero(M)
    ∂M⁻¹p∂q = zero(M)
    Ṁq̇ = @MVector zeros(T,size(M,2))
    ∂T∂qᵀ = @MVector zeros(T,size(M,2))
    (;S,x) = cf
    So = S(x(zeros(N)))
    Sg = ANCF.build_Sg(ancs,mass)
    Sps = [typeof(So)(S(x(loci[i].position))) for i in 1:num_of_loci]
    FlexibleBodyCoordinatesCache(
        pres_idx,
        free_idx,
        pres_idx_mask,
        free_idx_mask,
        cstr_idx,num_of_cstr,
        cf,e,ė,
        M,M⁻¹,
        ∂Mq̇∂q,∂M⁻¹p∂q,Ṁq̇,
        ∂T∂qᵀ,So,Sg,Sps)
end

struct FlexibleBodyState{N,M,T,cacheType} <: AbstractFlexibleBodyState{N,T}
    mass_locus_state::LocusState{N,M,T}
	"Positions of anchor points"
    loci_states::Vector{LocusState{N,M,T}}
	"Cache"
    cache::cacheType
end

function FlexibleBodyState(prop::FlexibleBodyProperty{N,T},
        ancs::ANCF.ANC,
        e::AbstractVector,
        ė=zero(e);
        pres_idx=Int[],
        cstr_idx=get_cstr_idx(ancs)) where {N,T}
    (;mass_locus,loci) = prop
    num_of_loci = length(loci)
    nc = ANCF.get_num_of_coords(ancs)
    cache = get_CoordinatesCache(
        prop,
        ancs,
        MVector{nc}(e),
        MVector{nc}(ė);
        pres_idx,cstr_idx
    )
    (;Sg,Sps) = cache
    mass_locus_state = LocusState(
        MVector{N}(Sg*e),
        MVector{N}(Sg*ė),
    )
    loci_states = [
        LocusState(
            MVector{N}(Sps[i]*e),
            MVector{N}(Sps[i]*ė)
        )
        for i in 1:num_of_loci
    ]
    FlexibleBodyState(
        mass_locus_state,
        loci_states,
        cache,
    )
end

struct FlexibleBody{N,M,T,cacheType,meshType} <: AbstractFlexibleBody{N,T}
	"Properties of a flexible body"
    prop::FlexibleBodyProperty{N,T}
	"State of a flexible body"
    state::FlexibleBodyState{N,M,T,cacheType}
    "Mesh"
    mesh::meshType
end

function FlexibleBody(prop,state)
    FlexibleBody(prop,state,nothing)
end


function body_state2coords_state(fb::FlexibleBody)
    (;e,ė) = fb.state.cache
    e,ė
end

function update_state!(
        state::FlexibleBodyState,
        cache::FlexibleBodyCoordinatesCache{<:ANCF.CoordinateFunctions},
        prop::AbstractBodyProperty,q,q̇)
    (;cache,mass_locus_state) = state
    (;e, ė, Sg) = cache
    e .= q
    ė .= q̇
    mul!(mass_locus_state.position, Sg, q)
    mul!(mass_locus_state.velocity, Sg, q̇)
end

function stretch_loci!(
    cache::FlexibleBodyCoordinatesCache,
	prop::AbstractBodyProperty,c)
end

function update_transformations!(
    cache::FlexibleBodyCoordinatesCache,
    state::FlexibleBodyState,
    prop::AbstractBodyProperty,q)
end

function update_loci_states!(
        state::FlexibleBodyState,
        cache::FlexibleBodyCoordinatesCache{<:ANCF.CoordinateFunctions},
        prop::AbstractBodyProperty,e,ė)
    (;loci_states) = state
    (;Sps) = cache
    for (i,locus_state) in enumerate(loci_states)
        mul!(locus_state.position, Sps[i], e)
        mul!(locus_state.velocity, Sps[i], ė)
    end
end

function generalize_force!(F,state::AbstractFlexibleBodyState)
	(;cache,mass_locus_state,loci_states) = state
	(;Sps,Sg,e,∂T∂eᵀ) = cache
	for (pid,locus_state) in enumerate(loci_states)
		mul!(F,transpose(Sps[pid]),locus_state.force,1,1)
	end
	mul!(F,transpose(Sg),mass_locus_state.force,1,1)
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
    (;mass_locus_state) = fb.state
    gravity_acceleration = get_gravity(fb)
    -transpose(mass_locus_state.position)*gravity_acceleration*mass
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
            mass_locus = SVector(subL/2,T(0),T(0))
            r̄p1 = SVector(T(0),T(0),T(0))
            r̄p2 = SVector(subL,T(0),T(0))
            loci = [
                r̄p1,r̄p2
            ]
            prop = FlexibleBodyProperty(
                i,
                :cable,
                mass,
                mass_locus,
                # length(loci),
                loci
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
    (;mass_locus_state,loci_states) = fb.state
    mass_locus_state.force .= 0
    foreach(loci_states) do locus_state
        locus_state.force .= 0
    end
end
