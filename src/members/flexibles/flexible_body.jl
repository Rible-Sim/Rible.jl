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

struct FlexibleBodyCoordinatesCache{fT,eT,inertia_cacheType,ArrayT}
    funcs::fT
    e::eT
    ė::eT
    inertia_cache::inertia_cacheType
    So::ArrayT
    Sg::ArrayT
    Sps::Vector{ArrayT}
end

function get_cstr_idx(ancs::ANCF.ANC)
	num_of_cstr = ANCF.get_num_of_cstr(ancs)
	collect(1:num_of_cstr)
end

function BodyCache(
        prop::FlexibleBodyProperty{N,T},
        ancs::ANCF.ANC,
        e::AbstractVector,
        ė=zero(e);
        pres_idx=Int[],
        cstr_idx=get_cstr_idx(ancs)
    ) where {N,T}
    (;mass,loci) = prop
    num_of_loci = length(loci)
    cf = ANCF.CoordinateFunctions(ancs,)
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
    inertia_cache = InertiaCache(
        M,
        M⁻¹,
        ∂Mq̇∂q,
        ∂M⁻¹p∂q,
        Ṁq̇,
        ∂T∂qᵀ,
    )
    FlexibleBodyCoordinatesCache(
        cf,e,ė,
        inertia_cache,
        So,Sg,Sps
    )
end

struct FlexibleBodyState{N,M,T} <: AbstractFlexibleBodyState{N,T}
    mass_locus_state::LocusState{N,M,T}
	"Positions of anchor points"
    loci_states::Vector{LocusState{N,M,T}}
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
    cache = BodyCache(
        prop,
        ancs,
        MVector{nc}(e),
        MVector{nc}(ė);
        pres_idx,cstr_idx
    )
    (;Sg,Sps) = cache
    mass_locus_state = LocusState(
        SVector{N}(Sg*e),
        SVector{N}(Sg*ė),
    )
    loci_states = [
        LocusState(
            SVector{N}(Sps[i]*e),
            SVector{N}(Sps[i]*ė)
        )
        for i in 1:num_of_loci
    ]
    coords = NonminimalCoordinates(
        ancs,pres_idx,cstr_idx
    )
    state = FlexibleBodyState(
        mass_locus_state,
        loci_states,
    )
    state, coords, cache
end

struct FlexibleBody{N,M,T,coordsType,cacheType,meshType} <: AbstractFlexibleBody{N,T}
	"Properties of a flexible body"
    prop::FlexibleBodyProperty{N,T}
	"State of a flexible body"
    state::FlexibleBodyState{N,M,T}
    "Coordinates Info"
    coords::NonminimalCoordinates{coordsType}
    "Cache"
    cache::cacheType
    "Mesh"
    mesh::meshType
end

function FlexibleBody(prop,state,coords,cache)
    FlexibleBody(prop,state,coords,cache,nothing)
end


function body_state2coords_state(fb::FlexibleBody)
    (;e,ė) = fb.cache
    e,ė
end

function lazy_update_state!(
        state::FlexibleBodyState,
        coords::NonminimalCoordinates{<:ANCF.ANC},
        cache::FlexibleBodyCoordinatesCache,
        prop::FlexibleBodyProperty,q,q̇)
    update_state!(state,coords,cache,prop,q,q̇)
end

function update_state!(
        state::FlexibleBodyState,
        coords::NonminimalCoordinates{<:ANCF.ANC},
        cache::FlexibleBodyCoordinatesCache,
        prop::FlexibleBodyProperty,q,q̇)
    (;mass_locus_state) = state
    (;e, ė, Sg) = cache
    e .= q
    ė .= q̇
    mass_locus_state.frame.position = Sg*q
    mass_locus_state.frame.velocity = Sg*q̇
end

function update_inertia_cache!(
        cache::FlexibleBodyCoordinatesCache,
        coords::NonminimalCoordinates{<:ANCF.ANC},
        prop::FlexibleBodyProperty,
        q,q̇
    )
    # constant mass matrix, no need to update
end


function stretch_loci!(
        coords::NonminimalCoordinates{<:ANCF.ANC},
        cache,
        prop::FlexibleBodyProperty,c
    )
    # to be implemented
end


function update_transformations!(
    coords::NonminimalCoordinates{<:ANCF.ANC},
    cache::FlexibleBodyCoordinatesCache,
    state::FlexibleBodyState,
    prop::FlexibleBodyProperty,q)
    # to be implemented
end

function update_loci_states!(
        state::FlexibleBodyState,
        coords::NonminimalCoordinates{<:ANCF.ANC},
        cache::FlexibleBodyCoordinatesCache,
        prop::FlexibleBodyProperty,e,ė)
    (;loci_states) = state
    (;Sps) = cache
    for (i,locus_state) in enumerate(loci_states)
        locus_state.frame.position = Sps[i]*e
        locus_state.frame.velocity = Sps[i]*ė
    end
end

function centrifugal_force!(F,state::AbstractFlexibleBodyState,cache)
	(;mass_locus_state,loci_states) = state
    (;∂T∂qᵀ) = cache.inertia_cache
	F .+= ∂T∂qᵀ
end

function mass_center_force!(F,state::AbstractFlexibleBodyState,cache)
	(;mass_locus_state,loci_states) = state
	(;Sps,Sg,e,) = cache
	mul!(F,transpose(Sg),mass_locus_state.force,1,1)
end

function concentrated_force!(F,state::AbstractFlexibleBodyState,cache)
	(;mass_locus_state,loci_states) = state
	(;Sps,Sg,e,) = cache
	for (pid,locus_state) in enumerate(loci_states)
		mul!(F,transpose(Sps[pid]),locus_state.force,1,1)
	end    
end

function strain!(F,state::AbstractFlexibleBodyState,cache)
	(;mass_locus_state,loci_states) = state
	(;Sps,Sg,e,) = cache
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
    -transpose(mass_locus_state.frame.position)*gravity_acceleration*mass
end

"""
Return the potential energy for strain of a flexible body
$(TYPEDSIGNATURES)
"""
function potential_energy_strain(fb::AbstractFlexibleBody)
    (;coords,cache) = fb
    (;e) = cache
    ancs = coords.nmcs
    ANCF.make_V(ancs)(e)
end

function subdivide(fb::FlexibleBody,nx=1,ny=1,nz=1)
    (;prop,state,coords,cache) = fb
    (;e,ė) = cache
    (;nmcs,pres_idx) = coords
    ancs = nmcs
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
            # cache = BodyCache(prop,ancs,𝐞)
            state, coords, cache = FlexibleBodyState(prop,subancs,esegs[i],ėsegs[i];pres_idx=sub_pres_idx)
            fb = FlexibleBody(prop,state, coords, cache)
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
