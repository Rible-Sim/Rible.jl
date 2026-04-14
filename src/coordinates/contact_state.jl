"""
ContactState Type.
$(TYPEDEF)
"""
struct MonoContactCoordinatesState{S <: AbstractCoordinatesState, CT} <: AbstractCoordinatesState
    state::S
    Λ::CT # normal contact impulses
    Γ::CT # tangential contact impulses
end

function Base.getproperty(st::MonoContactCoordinatesState, sym::Symbol)
    state = getfield(st, :state)
    if sym === :Λ
        return getfield(st, :Λ)
    elseif sym === :Γ
        return getfield(st, :Γ)
    else
        return getproperty(state, sym)
    end
end

function Base.setproperty!(st::MonoContactCoordinatesState, sym::Symbol, value)
    state = getfield(st, :state)
    if sym === :Λ
        return setfield!(st, :Λ, value)
    elseif sym === :Γ
        return setfield!(st, :Γ, value)
    else
        return setproperty!(state, sym, value)
    end
end

function Base.propertynames(st::MonoContactCoordinatesState)
    (propertynames(st.state)..., :Λ, :Γ)
end

"""
ContactState Type.
$(TYPEDEF)
"""
struct InnerContactCoordinatesState{S <: AbstractCoordinatesState, CT} <: AbstractCoordinatesState
    state::S
    Λ::CT # normal contact impulses
end

function Base.getproperty(st::InnerContactCoordinatesState, sym::Symbol)
    state = getfield(st, :state)
    if sym === :Λ
        return getfield(st, :Λ)
    else
        return getproperty(state, sym)
    end
end

function Base.setproperty!(st::InnerContactCoordinatesState, sym::Symbol, value)
    state = getfield(st, :state)
    if sym === :Λ
        return setfield!(st, :Λ, value)
    else
        return setproperty!(state, sym, value)
    end
end

function Base.propertynames(st::InnerContactCoordinatesState)
    (propertynames(st.state)..., :Λ)
end