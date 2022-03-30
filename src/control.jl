# include("PIDController.jl")
#
# using .PIDController
# import .PIDController: reset!

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
abstract type Heater <: AbstractController end
abstract type GroupedController <: AbstractController end
abstract type GroupedActuator <: GroupedController end
abstract type GroupedHeater <: GroupedController end

struct ManualActuator{CT,TT} <:AbstractController
    reg::CT
    traj::TT
end

struct ManualGangedActuators{CT,TT} <:GroupedActuator
    regs::CT
    traj::TT
end

struct ManualSerialActuators{CT,TT} <:GroupedActuator
    regs::CT
    traj::TT
end

struct ManualHeater{CT,HT,TT} <:Heater
    reg::CT
    heating_law::HT
    traj::TT
end

struct ManualSerialHeater{CT,HT,TT} <:GroupedHeater
    regs::CT
    heating_laws::HT
    trajs::TT
end

abstract type ControlScheme end

struct EmptyScheme <: ControlScheme end

# function record!(traj::ControlTrajectory,pid::PID)
#     @unpack ts,es,us = traj
#     push!(ts,pid.lastTime)
#     push!(es,pid.lastErr)
#     push!(us,pid.lastOutput)
# end

abstract type AbstractRegistor{T} end

struct SimpleRegistor{T} <: AbstractRegistor{T}
    id_string::Int
    original_value::T
end

struct SimpleRegistors{T} <: AbstractRegistor{T}
    id_cables::Vector{Int}
    original_values::Vector{T}
end

function ManualActuator(reg::SimpleRegistor{T}) where {T}
    ts = Vector{T}()
    es = Vector{T}()
    us = Vector{T}()
    traj = ControlTrajectory(ts,es,us)
    ManualActuator(reg,traj)
end

function ManualGangedActuators(regs::SimpleRegistors{T}) where {T}
    ts = Vector{T}()
    es = Vector{T}()
    us = Vector{T}()
    traj = ControlTrajectory(ts,es,us)
    ManualGangedActuators(regs,traj)
end

function ManualSerialActuators(regs::SimpleRegistors{T}) where {T}
    ts = Vector{Vector{T}}()
    es = Vector{Vector{T}}()
    us = Vector{Vector{T}}()
    trajs = ControlTrajectory(ts,es,us)
    ManualSerialActuators(regs,trajs)
end

function ManualHeater(reg::SimpleRegistor{T},heating_law) where {T}
    ts = Vector{T}()
    es = Vector{T}()
    us = Vector{T}()
    traj = ControlTrajectory(ts,es,us)
    ManualHeater(reg,heating_law,traj)
end

function ManualSerialHeater(regs::SimpleRegistors{T},heating_laws) where {T}
    ts = Vector{Vector{T}}()
    es = Vector{Vector{T}}()
    us = Vector{Vector{T}}()
    trajs = ControlTrajectory(ts,es,us)
    ManualSerialHeater(regs,heating_laws,trajs)
end

select_by_id(xs,id) = xs[findfirst((x)->x.id==id, xs)]

function reset!(ctrller::ManualActuator,tg)
    @unpack reg, traj = ctrller
    @unpack id_string, original_value = reg
    @unpack cables = tg
    s = select_by_id(cables,id_string)
    s.state.restlen = original_value
    reset!(traj)
end

function reset!(ctrller::GroupedActuator,tg)
    @unpack regs, trajs = ctrller
    @unpack id_cables, original_values = regs
    @unpack cables = tg
    for (id,original_value) in zip(id_cables,original_values)
        s = select_by_id(cables,id)
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



function actuate!(ctrller::ManualActuator,tg,u;inc=false)
    @unpack cables = tg
    @unpack reg, traj = ctrller
    @unpack id_string, original_value = reg
    @unpack us = traj
    s = select_by_id(cables,id_string)
    if inc
        s.state.restlen += u
    else
        s.state.restlen = original_value + u
    end
    # push!(us,u)
end

function actuate!(ctrller::ManualGangedActuators,tg,u;inc=false)
    @unpack cables = tg
    @unpack regs, traj = ctrller
    @unpack id_cables, original_values = regs
    @unpack us = traj
    s1 = select_by_id(cables,id_cables[1])
    s2 = select_by_id(cables,id_cables[2])
    if inc
        s1.state.restlen += u
        s2.state.restlen -= u
    else
        s1.state.restlen = original_values[1] + u
        s2.state.restlen = original_values[2] - u
    end
    # push!(us,u)
end

function actuate!(ctrller::ManualSerialActuators,tg,u;inc=false)
    @unpack cables = tg
    @unpack acts, trajs = ctrller
    @unpack id_cables, original_values = acts
    @unpack us = trajs
    for (id, original_value) in zip(id_cables,original_values)
        s = select_by_id(cables,id)
        if inc
            s.state.restlen += u
        else
            s.state.restlen = original_value + u
        end
    end
    # push!(us,u)
end

function heat!(ctrller::ManualHeater,tg,u;inc=false,abs=true)
    @unpack SMA_cables = tg.tensiles
    @unpack id_string, original_value, heating_law = ctrller.act
    s = select_by_id(SMA_cables,id_string)
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

@inline function heat!(ctrller::ManualSerialHeater,tg::TensegrityStructure,u;inc=false)
    heat!(ctrller,tg.tensiles.SMA_cables,u;inc)
end

function heat!(ctrller::ManualSerialHeater,
                SMA_cables::AbstractVector{<:SMACable},u;inc=false,abs=true)
    @unpack acts, heating_laws = ctrller
    @unpack id_cables, original_values = acts
    for (id, original_value, heating_law) in zip(id_cables,original_values,heating_laws)
        s = select_by_id(SMA_cables,id)
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

function heat!(s::SMACable,heating_law)
    s.law.F0, s.law.k = heating_law(s.state.temp)
end

function set_restlen!(tg,u)
    for (i,s) in enumerate(tg.cables)
        s.state.restlen = u[i]
    end
end
