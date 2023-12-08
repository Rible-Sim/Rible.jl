# include("PIDController.jl")
#
# using .PIDController
# import .PIDController: reset!
"""
所有？？超类。
"""
abstract type AbstractCoupler end
"""
所有？？制动超类。
"""
abstract type AbstractActuator{CT} end
struct Uncoupled <: AbstractCoupler end
struct Ganged <: AbstractCoupler end
struct Serial <: AbstractCoupler end

struct ManualActuator{CT<:AbstractCoupler,RT} <: AbstractActuator{CT}
    id::Int
    coupler::CT
    reg::RT
end

function ManualActuator(actid,id::Int,value::Number)
    ManualActuator(actid,Serial(),(ids=[id],values=[value]))
end

function ManualActuator(actid,ids::AbstractVector,values::AbstractVector,coupler=Serial())
    ManualActuator(actid,coupler,@eponymtuple(ids,values))
end

function actuate!(st,act::AbstractActuator{<:Uncoupled},μ::AbstractVector)
    (;cables) = st.tensiles
    (;reg) = act
    (;ids, values) = reg
    for (id, original_restlen) in zip(ids,values)
        cable = select_by_id(cables,id)
        cable.state.restlen = original_restlen + μ[id]
    end
end

function actuate!(st,act::AbstractActuator{<:Serial},μ::Number)
    (;cables) = st.tensiles
    (;reg) = act
    (;ids, values) = reg
    for (id, original_restlen) in zip(ids,values)
        cable = select_by_id(cables,id)
        cable.state.restlen = original_restlen + μ
    end
end

function actuate!(st,act::AbstractActuator{<:Ganged},μ::Number)
    (;cables) = st
    (;reg) = act
    (;ids, values) = reg
    cable1 = select_by_id(cables,ids[1])
    cable2 = select_by_id(cables,ids[2])
    cable1.state.restlen = values[1] + μ
    cable2.state.restlen = values[2] - μ
end

struct PrescribedActuator{MT,FT}
    id::Int
    manual::MT
    pres::FT
end

function actuate!(st,act::PrescribedActuator,t::Real)
    actuate!(st,act.manual,act.pres(t))
end
#
# struct ManualHeater{CT,HT,TT} <:Heater
#     reg::CT
#     heating_law::HT
#     traj::TT
# end
#
# struct ManualSerialHeater{CT,HT,TT} <:GroupedHeater
#     regs::CT
#     heating_laws::HT
#     trajs::TT
# end
#
#
# function ManualHeater(reg::SimpleRegistor{T},heating_law) where {T}
#     ts = Vector{T}()
#     restitution_coefficients = Vector{T}()
#     us = Vector{T}()
#     traj = ControlTrajectory(ts,restitution_coefficients,us)
#     ManualHeater(reg,heating_law,traj)
# end
#
# function ManualSerialHeater(regs::SimpleRegistors{T},heating_laws) where {T}
#     ts = Vector{Vector{T}}()
#     restitution_coefficients = Vector{Vector{T}}()
#     us = Vector{Vector{T}}()
#     trajs = ControlTrajectory(ts,restitution_coefficients,us)
#     ManualSerialHeater(regs,heating_laws,trajs)
# end

select_by_id(xs,id) = xs[findfirst((x)->x.id==id, xs)]

function actuate!(bot::Robot,friction_coefficients)
    (;st, hub) = bot
    (;actuators) = hub
    foreach(actuators) do actuator
        (;id) = actuator
        actuate!(st,actuator,friction_coefficients[id])
    end
end

# function heat!(tr::Robot,us;inc=false)
#     (;st, hub) = tr
#     (;heaters) = hub
#     for (heater,u) in zip(heaters,us)
#         heat!(heater,st,u;inc)
#     end
# end
#
# function heat!(ctrller::ManualHeater,st,u;inc=false,abs=true)
#     (;SMA_cables) = st.tensiles
#     (;id_string, original_value, heating_law) = ctrller.act
#     s = select_by_id(SMA_cables,id_string)
#     if abs
#         s.state.temp = u
#     else
#         if inc
#             s.state.temp += u
#         else
#             s.state.temp = original_value + u
#         end
#     end
#     heat!(s,heating_law)
# end
#
# @inline function heat!(ctrller::ManualSerialHeater,st::Structure,u;inc=false)
#     heat!(ctrller,st.tensiles.SMA_cables,u;inc)
# end
#
# function heat!(ctrller::ManualSerialHeater,
#                 SMA_cables::AbstractVector{<:SMADistanceSpringDamper},u;inc=false,abs=true)
#     (;acts, heating_laws) = ctrller
#     (;id_cables, original_values) = acts
#     for (id, original_value, heating_law) in zip(id_cables,original_values,heating_laws)
#         s = select_by_id(SMA_cables,id)
#         if abs
#             s.state.temp = u
#         else
#             if inc
#                 s.state.temp += u
#             else
#                 s.state.temp = original_value + u
#             end
#         end
#         heat!(s, heating_law)
#     end
# end
#
# function heat!(s::SMADistanceSpringDamper,heating_law)
#     s.law.F0, s.law.k = heating_law(s.state.temp)
# end

function set_restlen!(st,u)
    for (i,s) in enumerate(st.tensiles.cables)
        s.state.restlen = u[i]
    end
end
