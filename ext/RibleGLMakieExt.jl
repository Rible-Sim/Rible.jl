module RibleGLMakieExt

import Rible as RB
import Makie
import GLMakie as GM

function GM.draw_atomic(screen::GM.Screen, scene::Makie.Scene,vis::RB.Vis{Tuple{S}};) where {S<:RB.Apparatus}
end

end # module