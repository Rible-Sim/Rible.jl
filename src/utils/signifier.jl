
"""
Signifier — identity of a measurement point on a body.

Used by gauges/capta to index into `body.state.loci_states[pid]`.

$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
struct Signifier{bodyType}
    "Signifier of body"
    body::bodyType
    "Index of the locus"
    pid::Int
end

"""
Anchor — geometric attachment point for joints.

Holds the resolved position and axes for translation/rotation constraints,
independent of loci indexing.

$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
struct Anchor{bodyType,posType,trlAxesType,rotAxesType}
    "Body"
    body::bodyType
    "Position of the anchor in body-local frame"
    position::posType
    "Translational constraint axes"
    trl_axes::trlAxesType
    "Rotational constraint axes"
    rot_axes::rotAxesType
end

# ── Anchor constructors ──

"""
    Anchor(body, pid::Integer)

Construct from a single locus — position and axes all come from `body.prop.loci[pid]`.
"""
function Anchor(body, pid::Integer)
    locus = body.prop.loci[pid]
    Anchor(body, locus.position, locus.axes, locus.axes)
end

"""
    Anchor(body, position_pid::Integer, axes_pid::Integer)

Resolve position and axes from potentially different loci.
"""
function Anchor(body, position_pid::Integer, axes_pid::Integer)
    position_locus = body.prop.loci[position_pid]
    axes_locus = body.prop.loci[axes_pid]
    Anchor(body, position_locus.position, axes_locus.axes, axes_locus.axes)
end

"""
    Anchor(body, position, axes)

Direct geometry — same axes for translation and rotation.
"""
function Anchor(body, position, axes)
    Anchor(body, position, axes, axes)
end


"""
    Anchor(sig::Signifier)

Bridge from a Signifier — resolve geometry from the locus it identifies.
"""
function Anchor(sig::Signifier)
    Anchor(sig.body, sig.pid)
end

"""
Hen 2 Egg
$(TYPEDEF)
---
$(TYPEDFIELDS)
"""
struct Hen2Egg{henType,eggType}
    "hen/parent/predecessor"
    hen::henType
    "egg/child/successor"
    egg::eggType
end
