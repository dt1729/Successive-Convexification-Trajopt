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

function Linearisation(prevState::QuadrotorState, prevControl::QuadrotorControl)

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

function f(states::Matrix, control::Matrix)
    
    ans = zeros(6,100)
    ans[:,1] = states[:,1]
    
    for i in 2:100
    
        ans[1,i] = state[1, i-1] + control[1,i-1]*sin(state[4,i-1])*cos(state[5,i-1])*δ_t
        ans[2,i] = state[2, i-1] + control[1,i-1]*sin(state[4,i-1])*cos(state[5,i-1])*δ_t
        ans[3,i] = state[3, i-1] + control[1,i-1]*sin(state[5,i-1])*δ_t
        ans[4,i] = state[4, i-1] + control[2,i-1]*δ_t
        ans[5,i] = state[5, i-1] + control[3,i-1]*δ_t
        ans[6,i] = state[6, i-1] + 0.0
    
    end 

    return ans

end

function J(states::Matrix, control::Matrix, λ::Float64, D::Variable, W::Variable)
    tau = states - f(states, control)
    return sumsquares(P*evaluate(D)) + sumsquares(Q*evaluate(W)) + λ * maximum([sum(abs.(tau[:,i])) for i in 1:100])
end


function MPC_SuccessiveConvexification(referenceState::QuadrotorState, previousState::QuadrotorState, controlAction::QuadrotorControl) 
    
    # Issi ke andar integration bhi kar dena
    #        actual_change = last_nonlinear_cost - nonlinear_cost  # delta_J
    #        predicted_change = last_nonlinear_cost - linear_cost  # delta_L
    
    # Step 1 

    X_prev = zeros(6,100)
    U_prev = zeros(3,99)
    λ = 1e5
    Δ = 10.0
    δ_t = 0.1
    ρ_zero = 0.1 
    ρ_one  = 0.4
    ρ_two  = 0.5 
    α = 1.5
    count = 0
    
    while count < 20 

        D = Variable(6,100)
        W = Variable(3,99)
        v = Variable(6,100)

        
        A, B, D_t = Linearisation(prevState, prevControl)
        
        constraints = Constraint[X_prev[1,i] + D[1,i] >= 0  for i in 1:100]
        
        for i in 1:100

            for j in 1:3
                push!(constraints,0 <= X_prev[j,i] + D[j,i])
                push!(constraints,10 >= X_prev[j,i] + D[j,i])
            end

            if i < 100
                push!(constraints,-1 <=  U_prev[1,i] + W[1,i])
                push!(constraints, 1 >=  U_prev[1,i] + W[1,i])
                push!(constraints,-0.3 <= U_prev[2,i] + W[2,i])
                push!(constraints, 0.3 >= U_prev[2,i] + W[2,i])
                push!(constraints,-0.3 <= U_prev[3,i] + W[3,i])
                push!(constraints,0.3  >= U_prev[3,i] + W[3,i])
                push!(constraints, D[:,i+1] == A*D[:,i] + B*W[:,i] +  v[:,i] + D_t) 
                push!(constraints, norm(W[:,i], Inf) <= Δ)                   
            end

        end

        problem = minimize(sumsquares(P*D) + sumsquares(Q*W) + λ * sum(abs(v)), constraints)

        solve!(problem, SCS.Optimizer; silent_solver=True)
        
        ΔJ = J(X_prev + evaluate(D), U_prev + evaluate(W), λ, D, W) - J(X_prev, U_prev, λ, D, W)
        tt = evaluate(v)
        ΔL = J(X_prev + evaluate(D), U_prev + evaluate(W), λ, D, W) - (sumsquares(P*evaluate(D)) + sumsquares(Q*evaluate(W)) + λ * maximum([sum(abs.(tt[:,i])) for i in 1:100]))

        if ΔL <= 0.01 
            break 
        else
            ratio = ΔJ/ΔL
        end
        
        if ratio < ρ_zero
            Δ = Δ/α
            continue
        else
            X_prev = X_prev + evaluate(D)
            U_prev = U_prev + evaluate(W)
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

    outputControl = QuadrotorControl(U_prev[1,1], U_prev[2,1], U_prev[3,1])
    currentState  = QuadrotorState(0,0,0,0,0,0)

    currentState.x = state[1, 1] + U_prev[1,1]*sin(U_prev[2,1]*δ_t)*cos(U_prev[3,1]*δ_t)*δ_t
    currentState.y = state[2, 1] + U_prev[1,1]*sin(U_prev[2,1]*δ_t)*cos(U_prev[3,1]*δ_t)*δ_t
    currentState.z = state[3, 1] + U_prev[1,1]*sin(U_prev[3,1]*δ_t)*δ_t
    currentState.φ = state[4, 1] + U_prev[2,1]*δ_t
    currentState.θ = state[5, 1] + U_prev[3,1]*δ_t
    currentState.ψ = state[6, 1] + 0.0
        
    return currentState, outputControl
end

MPC_SuccessiveConvexification(QuadrotorState(1,1,1,0.1,0.1,0), QuadrotorState(0,0,0,0.0,0.0,0), QuadrotorControl(0,0,0))