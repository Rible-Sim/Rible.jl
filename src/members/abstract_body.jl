abstract type AbstractBody{N,T} end
abstract type AbstractBodyProperty{N,T} end
abstract type AbstractBodyState{N,T} end

update_body!(body::AbstractBody,q,q̇) = update_body!(body.state,body.state.cache,body.prop,q,q̇)
move_body!(body::AbstractBody,q,q̇)	= move_body!(body.state,body.state.cache,body.prop,q,q̇)
stretch_body!(body::AbstractBody,c) = stretch_body!(body.state.cache,body.prop,c)
update_transformations!(body::AbstractBody,q) = update_transformations!(body.state.cache,body.state,body.prop,q)
