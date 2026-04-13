# --- Extractors ---

"""
    extract(extractor, bot::Robot, )
    
Returns the state of interest `x`.
"""
function extract(extractor, bot::Robot) end


# --- Featurizers ---

"""
Return feature vector ϕ from state `x`.
$(TYPEDSIGNATURES)
"""
function features(featurizer, x) end

"""
Passes `x` through unchanged.
$(TYPEDEF)
"""
struct IdentityFeaturizer <: AbstractFeaturizer end

function features(::IdentityFeaturizer, x)
    return x
end


# --- Parametric Maps ---

"""
Return control action `u` from feature vector ϕ.
$(TYPEDSIGNATURES)
"""
function evaluate_paramfun end


# =============================================================================
# Gauge-Based Extractor (follows Objective pattern)
# =============================================================================

"""
Extracts state by measuring from gauges defined in `bot.hub.gauges`,
weighted by `gauges_weights`. Follows the same pattern as `Objective`.
$(TYPEDEF)
$(TYPEDFIELDS)
# Example
```julia
# Assuming bot.hub has 3 gauges defined
extractor = GaugesExtractor([1.0, 0.0, 1.0])  # use gauges 1 and 3

# All gauges with equal weight
extractor = GaugesExtractor(ones(3))
```

# Note
The gauges themselves are defined in `bot.hub.gauges` (same as used by Objective).
"""
struct GaugesExtractor{W<:AbstractVector} <: AbstractStateExtractor
    gauges_weights::W
end


"""
    extract(extractor, bot::Robot)

Extract state from robot at a specific instantaneous state.
"""
function extract(gex::GaugesExtractor, bot::Robot)
    (; structure, hub) = bot
    (; capta_gauges, coalition) = hub
    (; captum_gauge_id2sys_capta_idx) = coalition

    c = copy(hub.state.c)
    foreach(capta_gauges) do gauge
        c[captum_gauge_id2sys_capta_idx[gauge.id]] =  gex.gauges_weights[gauge.id] .* measure(structure, gauge)
    end
    return c
end


# =============================================================================
# Jacobian Support
# =============================================================================

"""
    extract_jacobian(extractor, bot::Robot)

Compute the Jacobian of the extracted state with respect to (q, q̇).
Returns (∂x/∂q, ∂x/∂q̇) where x is the extracted state.

Uses the same indexing pattern as `extract`: each gauge's Jacobian rows
are placed at indices given by `captum_gauge_id2sys_capta_idx[gauge.id]`.
"""
function extract_jacobian(gex::GaugesExtractor, bot::Robot)
    (;structure,hub) = bot
    (; capta_gauges, coalition) = hub
    (; captum_gauge_id2sys_capta_idx, num_of_capta) = coalition
    
    nq = get_num_of_full_coords(structure)
    ns = get_num_of_aux_var(structure)
    T = get_numbertype(structure)
    
    # Initialize Jacobian matrices with correct dimensions
    Jq = zeros(T, num_of_capta, nq)
    Jq̇ = zeros(T, num_of_capta, nq)
    Js = zeros(T, num_of_capta, ns)
    
    foreach(capta_gauges) do gauge
        w = gex.gauges_weights[gauge.id]
        idx = captum_gauge_id2sys_capta_idx[gauge.id]
        rows = length(idx)
        Jq_gauge = zeros(T, rows, nq)
        Jq̇_gauge = zeros(T, rows, nq)
        Js_gauge = zeros(T, rows, ns)
        measure_jacobian!(Jq_gauge, Jq̇_gauge, Js_gauge, structure, gauge.signifier, gauge.captum)
        Jq[idx, :] .= w .* Jq_gauge
        Jq̇[idx, :] .= w .* Jq̇_gauge
        Js[idx, :] .= w .* Js_gauge
    end
    
    return Jq, Jq̇, Js
end

"""
Policy that evaluates a time-dependent function `f(t)`.
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct TimeFunctionPolicy{F <: Function} <: AbstractTimePolicy
    f::F
end

# --- Generic Feedback Policy ---

"""
A modular feedback policy structure.
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct FeedbackPolicy{E <: AbstractStateExtractor, F <: AbstractFeaturizer, P <: AbstractParamFun} <: FeedbackDiscretePolicy
    extractor::E
    featurizer::F
    paramfun::P
end


"""
A discrete modular feedback policy structure.
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct DiscreteFeedbackPolicy{E <: AbstractStateExtractor, F <: AbstractFeaturizer, P <: AbstractParamFun} <: FeedbackDiscretePolicy
    extractor::E
    featurizer::F
    paramfun::P
end
