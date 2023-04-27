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

function TransformState(p1::QuadrotorState, tf::Transform3D{Float64})
    
    p1_coord = CartesianFrame3D("P1")
    
    newTf    = Transform3D(tt1, p1_coord, RotXYZ(p1.φ, p1.θ, p1.ψ), SVector(p1.x, p1.y, p1.z)) * tf
    
    angles   = RotXYZ(rotation(newTf))

    return QuadrotorState(newTf.mat[1, 4], newTf.mat[1, 5], newTf.mat[1,6], angles.theta1, angles.theta2, angles.theta3)
end

function distance(p1::QuadrotorState, p2::QuadrotorState)
    return sqrt((p1.x - p2.x)^2 + (p1.y - p2.y)^2 + (p1.z - p2.z)^2)
end

function TrajectoryTrackingLoop(referenceTrajectory::SVector, currentState::QuadrotorState)
    
    referenceState = QuadrotorState(0, 0, 0, 0, 0, 0)
    TransformedReferenceState = referenceTrajectory[0]
    count = 0
    base = CartesianFrame3D("base")
    Transforms = Vector{Transform3D{Float64}}()
    while !isCompleted
        
        currentState, prevControl = MPC_SuccessiveConvexification(TransformedReferenceState, QuadrotorState(0,0,0,0,0,0), prevControl)
        
        #Update local coordiante frame to the quadrotor position, save the position in global coordinate frame.
        
        tt1 = CartesianFrame3D("T"*string(count))

        if count == 0
      
            tempT = inv(Transform3D(tt1, base, RotXYZ(currentState.φ, currentState.θ, currentState.ψ), SVector(currentState.x, currentState.y, currentState.z)))
            push!(Transforms, Transform3D(tt1, base, RotXYZ(currentState.φ, currentState.θ, currentState.ψ), SVector(currentState.x, currentState.y, currentState.z)))

        else
      
            tempT = inv(Transforms[length(Transforms)])*inv(Transform3D(tt1, Transforms[length(Transforms)], RotXYZ(currentState.φ, currentState.θ, currentState.ψ), SVector(currentState.x, currentState.y, currentState.z)))
            push!(Transforms, Transform3D(tt1, Transforms[len(Transforms)], RotXYZ(currentState.φ, currentState.θ, currentState.ψ), SVector(currentState.x, currentState.y, currentState.z)))
   
        end

        
        if distance(TransformedreferenceTrajectory[i], currentState) <= thresholdDist
            i += 1 
        end
        
        TransformedReferenceState = TransformState(referenceTrajectory[i], Transforms[length(Trasnforms)], tt1)
        
        push!(quadrotorStateTrajectory, TransformState(currentState, TransformState, tt1))
        
        push!(quadrotorControlHistory, prevControl)
        
        count += 1
    end

end 

function MPC_SuccessiveConvexification(referenceState::QuadrotorState, previousState::QuadrotorState, controlAction::QuadrotorControl) 
    # Issi ke andar integration bhi kar dena


    outputControl = QuadrotorControl(0,0,0)

    return currentState, outputControl
end