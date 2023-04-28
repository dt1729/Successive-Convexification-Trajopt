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
#   i.) Write System equation A(t), B(t), D(t) - Done 
#   ii.) setup convex optimisation problem solver - Done
#   iii.) write trust region contraction code(L(k), J(k)) - Done
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

function Linearisation(time::Float64, prevState::QuadrotorState, prevControl::QuadrotorControl)

    A_t = zeros(Float64, 6, 6)
    B_t = zeros(Float64, 6, 3)
    D_t = zeros(Float64, 6, 1)

    A_t[1,3] = -prevControl[1]*sin(prevState[4])*cos(prevState[5])
    A_t[1,4] = -prevControl[1]*cos(prevState[4])*sin(prevState[5])
    A_t[2,3] = prevControl[1]*cos(prevState[4])*cos(prevState[5])
    A_t[2,4] = -prevControl[1]*sin(prevState[4])*sin(prevState[5])
    A_t[3,4] = prevControl[1]*cos(prevState[5])

    B_t[1,1] = cos(prevState[4])*cos(prevState[5])
    B_t[2,1] = sin(prevState[4])*cos(prevState[5])
    B_t[3,1] = sin(prevState[5])
    B_t[4,2] = 1
    B_t[5,3] = 1

    D_t[1,1] = -prevControl[1]*prevControl[2]*sin(prevState[4])*cos(prevState[5]) - prevControl[1]*prevControl[3]*cos(prevState[4])*cos(prevState[5])
    D_t[2,1] = prevControl[1]*prevControl[2]*cos(prevState[4])*cos(prevState[5]) - prevControl[1]*prevControl[3]*sin(prevState[4])*sin(prevState[5])
    D_t[3,1] = prevControl[1]*prevControl[3]*cos(prevState[5])

    return A_t, B_t, D_t
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
    #        actual_change = last_nonlinear_cost - nonlinear_cost  # delta_J
    #        predicted_change = last_nonlinear_cost - linear_cost  # delta_L
    
    # Step 1 
    X_prev = Variable(6,100)
    U_prev = Variable(3,99)
    λ = 1e5
    
    while count < 20 
        D = Variable(6,100)
        W = Variable(3,99)
        v = Variable(6,99)

        p = sumofsquares(P*D) + sumofsquares(Q*W) + λ * sum(abs(v))

        problem = minimize(p, constraints)
        
        constraints = Constraint[
           0 <=  X_prev[j,i] + D[j,i] <= 10.0 for j in 1:3 for i in 1:100,
          -1 <=  U_prev[1,i] + W[1,i] <= 1    for i in 1:99,
          -0.3 <= U_prev[2,i] + W[2,i] <= 0.3 for i in 1:99,
          -0.3 <= U_prev[3,i] + W[3,i] <= 0.3 for i in 1:99,
            D[:,i+1] == A*D[:,i] + B*W[:,i] +  v[:,i] + D for i in 1:99,
            max(w) <= Δ        
        ]
        
        solve!(problem, SCS.Optimizer; silent_solver=True)
        
        ΔJ = J(X_prev + D, U_prev + W) - J()
        ΔL = J(X_prev + D, U_prev + W) - L()

        if ΔL <= 0.01 
            ans = [X, U]
            break 
        else
            ratio = ΔJ/ΔL
        end
        
        if ratio < ρ_zero
            Δ = Δ/α
            continue
        else
            X_prev = X_prev + D
            U_prev = U_prev + W
            if ratio < ρ_one
                Δ = Δ/α
            else if ρ_1 <= ratio && ratio < ρ_two
                Δ = Δ
            else if ρ_2 <= ratio
                Δ = α * Δ
            end
        end
        Δ = max(Δ, Δ_min) 
        count  = count + 1
    end
    
    outputControl = QuadrotorControl(0,0,0)

    return currentState, outputControl
end