using StaticArrays
using Random
using LinearAlgebra
using Distances
using Infiltrator
using Plots
using CairoMakie
using JuMP, Ipopt, EAGO


############################################################################################
# TODO:
# 1.) Define trajectories and reference values.(3 should be enough)
# 2.) Write controller loop 
#   i.) Write position prediction retrieval code from the overall desired trajectory 
# 3.) Write MPC for the code
#   Successive Convexification
#   i.) Write System equation A(t), B(t), D(t) 
#   ii.) setup convex optimisation problem solver
#   iii.) write trust region contraction code(L(k), J(k))
############################################################################################

struct QuadrotorState
    x::Float64
    y::Float64
    z::Float64
    ψ::Float64
    θ::Float64
    φ::Float64
end

struct QuadrotorControl
    v::Float64
    ω::Float64
    α::Float64
end

function distance(p1::QuadrotorState, p2::QuadrotorState)
    return sqrt((p1.x - p2.x)^2 + (p1.y - p2.y)^2 + (p1.z - p2.z)^2)
end

function newRefState(referenceTrajectory::SVector, currentState::QuadrotorState)
    for i in eachindex(referenceTrajectory)
        if(distance(currentState, referenceTrajectory[i]) <  prev_dist)
            referenceState = referenceTrajectory[i]
            if(i != 0)
                return referenceState, i
            else
                return referenceState, 0
            end
        end
    end
end

function TrajectoryTrackingLoop(referenceTrajectory::SVector, currentState::QuadrotorState)
    referenceState = QuadrotorState(0,0,0,0,0,0)
    while !isCompleted
        currentState, prevControl = MPC_SuccessiveConvexification(referenceTrajectory[i], currentState, prevControl)
        if distance(referenceTrajectory[i], currentState) <= thresholdDist
            i += 1 
        end
        push!(quadrotorStateTrajectory, currentState)
        push!(quadrotorControlHistory, prevControl)
    end
end

function MPC_SuccessiveConvexification(referenceState::QuadrotorState, previousState::QuadrotorState, controlAction::QuadrotorControl)
    outputControl::QuadrotorControl
    return currentState, outputControl
end