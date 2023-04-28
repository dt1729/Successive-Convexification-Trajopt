using StaticArrays
using Random
using LinearAlgebra
# using Distances
using Infiltrator
# using Plots
# using CairoMakie
using Rotations, RigidBodyDynamics
using Convex, SCS, ECOS, Ipopt, Plots, JuMP
using DelimitedFiles
# using JuMP, Ipopt, EAGO

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
    φ::Float64
    θ::Float64
    ψ::Float64
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

    A_t[1,3] = -prevControl.v*sin(prevState.φ)*cos(prevState.θ)
    A_t[1,4] = -prevControl.v*cos(prevState.φ)*sin(prevState.θ)
    A_t[2,3] = prevControl.v*cos(prevState.φ)*cos(prevState.θ)
    A_t[2,4] = -prevControl.v*sin(prevState.φ)*sin(prevState.θ)
    A_t[3,4] = prevControl.v*cos(prevState.θ)

    B_t[1,1] = cos(prevState.φ)*cos(prevState.θ)
    B_t[2,1] = sin(prevState.φ)*cos(prevState.θ)
    B_t[3,1] = sin(prevState.θ)
    B_t[4,2] = 1
    B_t[5,3] = 1

    D_t[1,1] = -prevControl.v*prevControl.ω*sin(prevState.φ)*cos(prevState.θ) - prevControl.v*prevControl.α*cos(prevState.φ)*cos(prevState.θ)
    D_t[2,1] = prevControl.v*prevControl.ω*cos(prevState.φ)*cos(prevState.θ) - prevControl.v*prevControl.α*sin(prevState.φ)*sin(prevState.θ)
    D_t[3,1] = prevControl.v*prevControl.α*cos(prevState.θ)

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

function f(state::Matrix, control::Matrix)
    δ_t = 0.1
    ans = zeros(6,100)
    ans[:,1] = state[:,1]
    
    for i in 2:100
    
        ans[1,i] = state[1, i-1] + control[1,i-1]*sin(state[4,i-1])*cos(state[5,i-1])*δ_t
        ans[2,i] = state[2, i-1] + control[1,i-1]*sin(state[4,i-1])*cos(state[5,i-1])*δ_t
        ans[3,i] = state[3, i-1] + control[1,i-1]*sin(state[5,i-1])*δ_t
        ans[4,i] = state[4, i-1] + control[2,i-1]*δ_t
        ans[5,i] = state[5, i-1] + control[3,i-1]*δ_t
        ans[6,i] = 0.0
    
    end 

    return ans

end

function J(states::Matrix, control::Matrix, λ::Float64, D::Variable, W::Variable)
    tau = states - f(states, control)
    v1 = 0
    v2 = 0
    for i in 1:100
        v1 += D.value[:,i]'*P*D.value[:,i]
        if i <100 
            v2 += W.value[:,i]'*Q*W.value[:,i]
        end 
    end
    v3 =  λ * maximum([sum(abs.(tau[:,i])) for i in 1:100])
    return v1 + v2 + v3
end

function VanillaOptimisation(referenceState::QuadrotorState, previousState::QuadrotorState, controlAction::QuadrotorControl)
    X_prev = zeros(6,100)
    U_prev = zeros(3,99)
    λ = 1e5
    Δ = 10.0
    Δ_min = 0.1
    δ_t = 0.1
    ρ_zero = 0.1 
    ρ_1  = 0.25
    ρ_2  = 0.9
    α = 2.5
    count = 0

    P = 1.0*I(6)
    Q = 3.0*I(3)
    
    D = Variable(6,100)
    W = Variable(3,99)
    v = Variable(6,100)

    
    A, B, D_t = Linearisation(previousState, controlAction)
    
    constraints = Constraint[0 <= D[1,i] for i in 1:100]
    for i in 1:100
        push!(constraints,10 >= D[1,i])
    end
    
    for i in 1:100

        for j in 2:3
            push!(constraints,0 <=  D[j,i])
            push!(constraints,10 >= D[j,i])
        end

        if i < 100
            push!(constraints, -10 <=  W[1,i])
            push!(constraints, 1 >=   W[1,i])
            push!(constraints,-0.15 <= W[2,i])
            push!(constraints, 0.15 >= W[2,i])
            push!(constraints,-0.15 <= W[3,i])
            push!(constraints,0.15 >= W[3,i])
            push!(constraints, D[:,i+1] - D[:,i] == (A*D[:,i] + B*W[:,i])*δ_t) 
        end

    end

    #OBSTACLE CONSTRAINTS: 
    # for i in 1:100
    #     push!(constraints,5 <= (X_prev[1,i] + D[1,i]))
    #     push!(constraints,5 <= (X_prev[2,i] + D[2,i]))
    # end

    ref_state_vec = zeros(Float64, 6,1)
    ref_state_vec[1,1] = referenceState.x
    ref_state_vec[2,1] = referenceState.y
    ref_state_vec[3,1] = referenceState.z
    ref_state_vec[4,1] = referenceState.φ
    ref_state_vec[5,1] = referenceState.θ
    ref_state_vec[6,1] = referenceState.ψ
    push!(constraints, D[:,100] == ref_state_vec)

    problem = minimize(sumsquares(P*D) + sumsquares(Q*W))

    solve!(problem, Ipopt.Optimizer; silent_solver=true)
    println("Final pt: ", D.value)
end

function TrajectoryOptimisation(referenceState::QuadrotorState, previousState::QuadrotorState, controlAction::QuadrotorControl) 
    
    # Step 1 

    X_prev = zeros(6,100)
    U_prev = zeros(3,99)
    λ = 1e5
    Δ = 10.0
    Δ_min = 0.1
    δ_t = 0.1
    ρ_zero = 0.1 
    ρ_1  = 0.25
    ρ_2  = 0.9
    α = 2.5
    count = 0

    P = 1.0*I(6)
    Q = 3.0*I(3)
    
    while count < 5 

        D = Variable(6,100)
        W = Variable(3,99)
        v = Variable(6,100)

        
        A, B, D_t = Linearisation(previousState, controlAction)
        constraints = Constraint[0 <= X_prev[1,i] + D[1,i] for i in 1:100]
        for i in 1:100
            push!(constraints,10 >= X_prev[1,i] + D[1,i])
        end
        
        for i in 1:100

            for j in 2:3
                push!(constraints,0 <= X_prev[j,i] + D[j,i])
                push!(constraints,10 >= X_prev[j,i] + D[j,i])
            end

            if i < 100
                push!(constraints, -1 <=  U_prev[1,i] + W[1,i])
                push!(constraints, 1 >=  U_prev[1,i] + W[1,i])
                push!(constraints,-0.3 <= U_prev[2,i] + W[2,i])
                push!(constraints, 0.3 >= U_prev[2,i] + W[2,i])
                push!(constraints,-0.3 <= U_prev[3,i] + W[3,i])
                push!(constraints,0.3  >= U_prev[3,i] + W[3,i])
                push!(constraints, D[:,i+1] - D[:,i] == (A*D[:,i] + B*W[:,i] +  v[:,i] + D_t)*δ_t) 
                push!(constraints, norm(W[:,i], Inf) <= Δ)
                push!(constraints, norm(v[:,i], Inf) <= Δ)
            end

        end

        #OBSTACLE CONSTRAINTS: 
        # for i in 1:100
        #     push!(constraints,5 <= (X_prev[1,i] + D[1,i]))
        #     push!(constraints,5 <= (X_prev[2,i] + D[2,i]))
        # end

        ref_state_vec = zeros(Float64, 6,1)
        ini_state_vec = zeros(Float64, 6,1)
        ref_state_vec[1,1] = referenceState.x
        ref_state_vec[2,1] = referenceState.y
        ref_state_vec[3,1] = referenceState.z
        ref_state_vec[4,1] = referenceState.φ
        ref_state_vec[5,1] = referenceState.θ
        ref_state_vec[6,1] = referenceState.ψ

        ini_state_vec[1,1] = previousState.x
        ini_state_vec[2,1] = previousState.y
        ini_state_vec[3,1] = previousState.z
        ini_state_vec[4,1] = previousState.φ
        ini_state_vec[5,1] = previousState.θ
        ini_state_vec[6,1] = previousState.ψ
        push!(constraints, X_prev[:,100] + D[:,100] == ref_state_vec)
        # push!(constraints, X_prev[1,1] + D[1,1] <= ini_state_vec[1,1] + 0.5)
        # push!(constraints, X_prev[2,1] + D[2,1] <= ini_state_vec[2,1] + 0.5)
        # push!(constraints, X_prev[3,1] + D[3,1] <= ini_state_vec[3,1] + 0.5)
        problem = minimize(sumsquares(P*D) + sumsquares(Q*W) + λ * opnorm(v,1), constraints)

        @time solve!(problem, ECOS.Optimizer; silent_solver=false)
        
        ΔJ = J(X_prev + D.value, U_prev + W.value, λ, D, W) - J(X_prev, U_prev, λ, D, W)
        tt = v.value

        v1 = 0
        v2 = 0

        for i in 1:100
            v1 += D.value[:,i]'*P*D.value[:,i]
            if i <100 
                v2 += W.value[:,i]'*Q*W.value[:,i]
            end 
        end
            
        ΔL = J(X_prev + D.value, U_prev + W.value, λ, D, W) - (v1 + v2 + λ * maximum([sum(abs.(tt[:,i])) for i in 1:100]))
        println("Delta: ", Δ)
        if ΔL <= 0.01
            Δ = max(Δ, Δ_min) 
            
            count  = count + 1    
            continue 
        else
            ratio = ΔJ/ΔL
        end
        
        if ratio < ρ_zero
            Δ = Δ/α
            Δ = max(Δ, Δ_min) 
            println("Here")
            count  = count + 1    
            continue
        else
            X_prev = X_prev + D.value
            U_prev = U_prev + W.value
            if ratio < ρ_1
                Δ = Δ/α
            elseif ρ_1 <= ratio && ratio < ρ_2
                Δ = Δ
            elseif ρ_2 <= ratio
                Δ = α * Δ
            end
        end
        Δ = max(Δ, Δ_min) 
        println("Here")
        count = count + 1

    end

    outputControl = QuadrotorControl(U_prev[1,1], U_prev[2,1], U_prev[3,1])
    println("Controls: ", X_prev[1,100]," ", X_prev[2,100]," ", X_prev[3,100], " ", X_prev[4,100]," ", X_prev[5,100]," ", X_prev[6,100])
    writedlm("States_landing.txt", X_prev', ",")
    writedlm("Controls_landing.txt", U_prev', ",")
    curr_x = previousState.x + U_prev[1,1]*sin(U_prev[2,1]*δ_t)*cos(U_prev[3,1]*δ_t)*δ_t
    curr_y = previousState.y + U_prev[1,1]*sin(U_prev[2,1]*δ_t)*cos(U_prev[3,1]*δ_t)*δ_t
    curr_z = previousState.z + U_prev[1,1]*sin(U_prev[3,1]*δ_t)*δ_t
    curr_φ = previousState.φ + U_prev[2,1]*δ_t
    curr_θ = previousState.θ + U_prev[3,1]*δ_t
    curr_ψ = previousState.ψ + 0.0
    currentState  = QuadrotorState(curr_x,curr_y,curr_z,curr_φ,curr_θ,curr_ψ)
        
    return currentState, outputControl
end

referenceState =  QuadrotorState(3,4,5,0.1,0.1,0.0)
previousState  =  QuadrotorState(0,0,0,0.0,0.0,0)
prevControl    =  QuadrotorControl(0.01,0.001,0.001)
P = 1.0*I(6)
Q = 3.0*I(3)

# VanillaOptimisation(referenceState, previousState, prevControl)
previousState, prevControl = TrajectoryOptimisation(referenceState, previousState, prevControl)