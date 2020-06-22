include("PIDController.jl")

using .PIDController
import .PIDController: reset!

struct ControlTrajectory{tType,eType,uType}
    ts::Vector{tType}
    es::Vector{eType}
    us::Vector{uType}
end

function ControlTrajectory()
    ts = Vector{Float64}()
    es = Vector{Float64}()
    us = Vector{Float64}()
    ControlTrajectory(ts,es,us)
end

function reset!(traj::ControlTrajectory)
    @unpack ts,es,us = traj
    resize!(ts,0)
    resize!(es,0)
    resize!(us,0)
end

struct ControlHub{ctrlType,trajType}
    ctrls::ctrlType
    trajs::trajType
end

function ControlHub(ctrls)
    trajs = [ControlTrajectory() for ctrid in eachindex(ctrls)]
    ControlHub(ctrls,trajs)
end

function record!(traj::ControlTrajectory,pid::PID)
    @unpack ts,es,us = traj
    push!(ts,pid.lastTime)
    push!(es,pid.lastErr)
    push!(us,pid.lastOutput)
end

function reset!(hub::ControlHub)
    @unpack ctrls,trajs = hub
    reset!.(ctrls)
    reset!.(trajs)
end

struct Actuator{N,T}
    strings::Vector{SString{N,T}}
end

function reset!(actuator::Actuator)
    @unpack strings = actuator
    str1 = strings[1]
    str2 = strings[2]
    str1.state.restlen = str1.original_restlen
    str2.state.restlen = str2.original_restlen
end

function actuate!(tgstruct::TensegrityStructure,us;inc=false)
    for (actuator,u) in zip(tgstruct.actuators,us)
        actuate!(actuator,u,inc=inc)
    end
end

function actuate!(actuator::Actuator,u;inc=false)
    @unpack strings = actuator
    str1 = strings[1]
    str2 = strings[2]
    if inc
        str1.state.restlen += u
        str2.state.restlen -= u
    else
        str1.state.restlen = str1.original_restlen + u
        str2.state.restlen = str2.original_restlen - u
    end
end
