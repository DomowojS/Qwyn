using Optimization
using OptimizationBBO
using Random
using LinearAlgebra

# Mutable struct to hold optimization parameters
mutable struct OptimizationParams
    a::Vector{Float64}
    b::Vector{Float64}
    c::Vector{Float64}
    d::Float64
    e::Float64
end

# Objective function with correct signature
function objective_function(x, params::OptimizationParams)
    # Create a copy of params to modify
    local_params = deepcopy(params)
    
    # Modify vectors a and b using design variable x
    local_params.a .= x .* local_params.a
    local_params.b .= x .+ local_params.b

    # Compute arbitrary scalar objective
    return sum(local_params.a.^2) + dot(local_params.b, local_params.c) + local_params.d * norm(x) + local_params.e
end

# Constraint function (inequality constraint)
function cons!(g, x, p::OptimizationParams)
    # Example constraint: sum of x elements should be less than 5
    g[1] = sum(x) - 5.0
    return g
end

# Main optimization setup
function run_optimization()
    # Initialize random seed for reproducibility
    Random.seed!(42)

    # Initial random design variable
    x0 = rand(3)

    # Parameter struct initialization
    params = OptimizationParams(
        [1.0, 2.0, 3.0],   # a
        [4.0, 5.0, 6.0],   # b
        [7.0, 8.0, 9.0],   # c
        2.5,               # d
        1.5                # e
    )

    # Define lower and upper bounds for x
    lower_bounds = [-2.0, -2.0, -2.0]
    upper_bounds = [2.0, 2.0, 2.0]

    # Create optimization function
    opt_func = OptimizationFunction(
        objective_function, 
        Optimization.AutoZygote(),
        #cons = cons!
    )

    # Problem definition
    prob = OptimizationProblem(opt_func, x0, params, 
        lb = lower_bounds, 
        ub = upper_bounds,
        #lcons = [0.0],
        #ucons = [Inf]
    )

    # Solve using BBO algorithm
    sol = solve(prob, BBO(), 
        maxiters = 100, 
        abstol = 1e-8, 
        reltol = 1e-8
    )

    # Print results
    println("Optimal solution: ", sol.u)
    println("Optimal objective value: ", sol.objective)
end

# Run the optimization
run_optimization()