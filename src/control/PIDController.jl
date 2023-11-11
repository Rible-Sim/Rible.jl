module PIDController
export PID,update!,tune!,reset!
using Parameters
mutable struct PID{T,limitType}
    Kp::T
    Ki::T
    Kd::T
    setpoint::T
    dt::T
    output_limits::limitType
    #
    input::T
    output::T
    pTerm::T
    iTerm::T
    dTerm::T
    lastErr::T
    lastInput::T
    lastOutput::T
    lastTime::T
    currentTime::T
end
function PID(Kp,Ki,Kd;setpoint,dt,initialInput=zero(Kp))
    @assert typeof(Kp)==typeof(Ki)==typeof(Kd)==typeof(setpoint)
    T = typeof(Kp)
    output_limits = nothing
    #
    input = zero(T)
    output = zero(T)
    pTerm = zero(T)
    iTerm = zero(T)
    dTerm = zero(T)
    lastErr = zero(T)
    lastInput = initialInput
    lastOutput = zero(T)
    lastTime = -dt
    currentTime = zero(T)
    PID(Kp,Ki,Kd,setpoint,
                dt,output_limits,
                input,output,
                pTerm,iTerm,dTerm,
                lastErr,lastInput,lastOutput,
                lastTime,
                currentTime)
end

function (pid::PID)(arg...)
    update!(pid,arg...)
end

function reset!(pid::PID)
    pid.input = 0.0
    pid.output = 0.0
    pid.pTerm = 0.0
    pid.iTerm = 0.0
    pid.dTerm = 0.0
    pid.lastErr = 0.0
    pid.lastInput = 0.0
    pid.lastOutput = 0.0
    pid.lastTime = -pid.dt
    pid.currentTime = 0.0
end
# function update!(pid,input,t)
#     pid.input = input
#     pid.currentTime = t
#     pid.dt = t - pid.lastTime
#     @assert pid.dt > 0.0
#     update!(pid)
# end

function update!(pid,input,t)
    if t - pid.lastTime >= pid.dt
        pid.input = input
        pid.currentTime = t
        update!(pid)
    else
        return pid.lastOutput
    end
end
#
# function update!(pid,input;dt)
#     @assert dt > 0.0
#     pid.input = input
#     pid.dt = dt
#     pid.currentTime = pid.lastTime + pid.dt
#     update!(pid)
# end

function update!(pid,input)
    pid.input = input
    pid.currentTime = pid.lastTime + pid.dt
    update!(pid)
end

function update!(pid)
    (;Kp,Ki,Kd,setpoint) = pid
    (;input,dt,output_limits) = pid
    # Compute all the working error variables
    err = setpoint - input
    # Compute PID every term
    dInput = input-pid.lastInput
    pid.pTerm = Kp*err
    pid.iTerm += Ki*err*dt
    pid.dTerm = -Kd*dInput/dt # Because setpoint doesn't change
    # Compute PID Output
    pid.output = pid.pTerm + pid.iTerm + pid.dTerm
    if typeof(output_limits)!=Nothing
        pid.output = clamp(pid.output,output_limits[1],output_limits[2])
    end
    # Remember some variables for next time
    pid.lastErr = err
    pid.lastOutput = pid.output
    pid.lastInput = pid.input
    pid.lastTime = pid.currentTime
    return pid.output
end
# function update!(pid,::Val{:OnMeasurement})
#     (;Kp,Ki,Kd,setpoint) = pid
#     (;input,dt,output_limits) = pid
#     # How long since we last calculated
#     @assert dt > 0.0
#     # Compute all the working error variables
#     err = setpoint - input
#     dInput = input-pid.lastInput
#     pid.pTerm -= Kp*dInput
#     pid.iTerm += Ki*err*dt
#     pid.dTerm = -Kd*dInput/dt
#     # Compute PID Output
#     pid.output = pid.pTerm + pid.iTerm + pid.dTerm
#     if typeof(output_limits)!=Nothing
#         pid.output = clamp(pid.output,output_limits[1],output_limits[2])
#     end
#     # Remember some variables for next time
#     pid.lastOutput = pid.output
#     pid.lastInput = pid.input
#     pid.lastTime = pid.currentTime
#     return pid.output
# end

function tune!(pid::PID,p,i,d)
    reset!(pid)
    pid.Kp = p
    pid.Ki = i
    pid.Kd = d
end

function tune!(pid::PID;p=pid.Kp,i=pid.Ki,d=pid.Kd)
    reset!(pid)
    pid.Kp = p
    pid.Ki = i
    pid.Kd = d
end
end
