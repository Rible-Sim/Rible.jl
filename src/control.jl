include("PIDController.jl")

using .PIDController
import .PIDController: reset!

struct ControlTrajectory{tType,eType,uType}
    ts::Vector{tType}
    es::Vector{eType}
    us::Vector{uType}
end

function reset!(traj::ControlTrajectory)
    @unpack ts,es,us = traj
    resize!(ts,0)
    resize!(es,0)
    resize!(us,0)
end

abstract type AbstractController end
abstract type Heating <: AbstractController end
abstract type GroupedController <: AbstractController end
abstract type GroupedActuation <: GroupedController end
abstract type GroupedHeating <: GroupedController end

struct ManualActuation{CT,TT} <:AbstractController
    act::CT
    traj::TT
end

struct ManualGangedActuation{CT,TT} <:GroupedActuation
    acts::CT
    traj::TT
end

struct ManualSerialActuation{CT,TT} <:GroupedActuation
    acts::CT
    traj::TT
end

struct ManualHeating{CT,HT,TT} <:Heating
    act::CT
    heating_law::HT
    traj::TT
end

struct ManualSerialHeating{CT,HT,TT} <:GroupedHeating
    acts::CT
    heating_laws::HT
    trajs::TT
end

abstract type ControlScheme end

struct EmptyScheme <: ControlScheme end

function record!(traj::ControlTrajectory,pid::PID)
    @unpack ts,es,us = traj
    push!(ts,pid.lastTime)
    push!(es,pid.lastErr)
    push!(us,pid.lastOutput)
end

abstract type AbstractActuator{T} end

struct SimpleActuator{T} <: AbstractActuator{T}
    id_string::Int
    original_value::T
end

struct SimpleActuators{T} <: AbstractActuator{T}
    id_strings::Vector{Int}
    original_values::Vector{T}
end

function ManualActuation(act::SimpleActuator{T}) where {T}
    ts = Vector{T}()
    es = Vector{T}()
    us = Vector{T}()
    traj = ControlTrajectory(ts,es,us)
    ManualActuation(act,traj)
end

function ManualGangedActuation(acts::SimpleActuators{T}) where {T}
    ts = Vector{T}()
    es = Vector{T}()
    us = Vector{T}()
    traj = ControlTrajectory(ts,es,us)
    ManualGangedActuation(acts,traj)
end

function ManualSerialActuation(acts::SimpleActuators{T}) where {T}
    ts = Vector{Vector{T}}()
    es = Vector{Vector{T}}()
    us = Vector{Vector{T}}()
    trajs = ControlTrajectory(ts,es,us)
    ManualSerialActuation(acts,trajs)
end

function ManualHeating(act::SimpleActuator{T},heating_law) where {T}
    ts = Vector{T}()
    es = Vector{T}()
    us = Vector{T}()
    traj = ControlTrajectory(ts,es,us)
    ManualHeating(act,heating_law,traj)
end

function ManualSerialHeating(acts::SimpleActuators{T},heating_laws) where {T}
    ts = Vector{Vector{T}}()
    es = Vector{Vector{T}}()
    us = Vector{Vector{T}}()
    trajs = ControlTrajectory(ts,es,us)
    ManualSerialHeating(acts,heating_laws,trajs)
end

select_by_id(xs,id) = xs[findfirst((x)->x.id==id, xs)]

function reset!(ctrller::ManualActuation,tg)
    @unpack act, traj = ctrller
    @unpack id_string, original_value = act
    @unpack strings = tg
    s = select_by_id(strings,id_string)
    s.state.restlen = original_value
    reset!(traj)
end

function reset!(ctrller::GroupedActuation,tg)
    @unpack acts, trajs = ctrller
    @unpack id_strings, original_values = acts
    @unpack strings = tg
    for (id,original_value) in zip(id_strings,original_values)
        s = select_by_id(strings,id)
        s.state.restlen = original_value
    end
    reset!(traj)
end

function actuate!(tr::TensegrityRobot,us;inc=false)
    @unpack tg, hub = tr
    @unpack actuators = hub
    for (actuator,u) in zip(actuators,us)
        actuate!(actuator,tg,u;inc)
    end
end

function heat!(tr::TensegrityRobot,us;inc=false)
    @unpack tg, hub = tr
    @unpack heaters = hub
    for (heater,u) in zip(heaters,us)
        heat!(heater,tg,u;inc)
    end
end



function actuate!(ctrller::ManualActuation,tg,u;inc=false)
    @unpack strings = tg
    @unpack act, traj = ctrller
    @unpack id_string, original_value = act
    @unpack us = traj
    s = select_by_id(strings,id_string)
    if inc
        s.state.restlen += u
    else
        s.state.restlen = original_value + u
    end
    # push!(us,u)
end

function actuate!(ctrller::ManualGangedActuation,tg,u;inc=false)
    @unpack strings = tg
    @unpack acts, traj = ctrller
    @unpack id_strings, original_values = acts
    @unpack us = traj
    s1 = select_by_id(strings,id_strings[1])
    s2 = select_by_id(strings,id_strings[2])
    if inc
        s1.state.restlen += u
        s2.state.restlen -= u
    else
        s1.state.restlen = original_values[1] + u
        s2.state.restlen = original_values[2] - u
    end
    # push!(us,u)
end

function actuate!(ctrller::ManualSerialActuation,tg,u;inc=false)
    @unpack strings = tg
    @unpack acts, trajs = ctrller
    @unpack id_strings, original_values = acts
    @unpack us = trajs
    for (id, original_value) in zip(id_strings,original_values)
        s = select_by_id(strings,id)
        if inc
            s.state.restlen += u
        else
            s.state.restlen = original_value + u
        end
    end
    # push!(us,u)
end

function heat!(ctrller::ManualHeating,tg,u;inc=false,abs=true)
    @unpack SMA_strings = tg.tensiles
    @unpack id_string, original_value, heating_law = ctrller.act
    s = select_by_id(SMA_strings,id_string)
    if abs
        s.state.temp = u
    else
        if inc
            s.state.temp += u
        else
            s.state.temp = original_value + u
        end
    end
    heat!(s,heating_law)
end

@inline function heat!(ctrller::ManualSerialHeating,tg::TensegrityStructure,u;inc=false)
    heat!(ctrller,tg.tensiles.SMA_strings,u;inc)
end

function heat!(ctrller::ManualSerialHeating,
                SMA_strings::AbstractVector{<:SMAString},u;inc=false,abs=true)
    @unpack acts, heating_laws = ctrller
    @unpack id_strings, original_values = acts
    for (id, original_value, heating_law) in zip(id_strings,original_values,heating_laws)
        s = select_by_id(SMA_strings,id)
        if abs
            s.state.temp = u
        else
            if inc
                s.state.temp += u
            else
                s.state.temp = original_value + u
            end
        end
        heat!(s, heating_law)
    end
end

function heat!(s::SMAString,heating_law)
    s.law.F0, s.law.k = heating_law(s.state.temp)
end

function set_restlen!(tgstruct,u)
    for (i,s) in enumerate(tgstruct.strings)
        s.state.restlen = u[i]
    end
end
