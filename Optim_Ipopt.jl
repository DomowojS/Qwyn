using OptimizationMOI
using Ipopt
using Random
using LinearAlgebra

# Define a mutable struct to hold optimization parameters
mutable struct OptimizationParams
    a::Vector{Float64}
    b::Vector{Float64}
    c::Vector{Float64}
    d::Float64
    e::Float64
end

# Objective function
function objective_func!(params::OptimizationParams, x)
    # Modify vectors a and b based on x
    params.a .= x[1:3] .* 2.0
    params.b .= x[4:6] .+ 1.0

    # Compute an arbitrary equation using all struct variables
    result = sum(params.a .* params.b) + 
             dot(params.c, params.a) + 
             params.d * maximum(params.b) - 
             params.e * minimum(params.a)
    
    return result
end

# Constraint function
function constraint!(g, x, p)
    # Inequality constraint: sum of x elements should be less than 4
    g[1] = sum(x) - 4.0
end

# Main optimization function
function run_optimization()
    # Set random seed for reproducibility
    Random.seed!(42)

    # Initialize optimization parameters
    params = OptimizationParams(
        zeros(3),   # a
        zeros(3),   # b
        [1.5, 2.5, 3.5],  # c
        2.0,        # d
        1.0         # e
    )

    # Initial design variable (random values)
    x0 = rand(6)

    # Define lower and upper bounds
    # Bounds for x: between -5 and 5, with first 3 elements more constrained
    lb = [-2.0, -2.0, -2.0, -5.0, -5.0, -5.0]
    ub = [2.0, 2.0, 2.0, 5.0, 5.0, 5.0]

    # Define the optimization function
    opt_func = OptimizationFunction(objective_func!, 
        cons! = constraint!, 
        num_cons = 1
    )

    # Set up the optimization problem
    prob = OptimizationProblem(
        opt_func, 
        x0, 
        params;
        lb = lb, 
        ub = ub
    )

    # Solve using Ipopt
    solver = Ipopt.Optimizer()
    sol = solve(prob, solver)

    # Print results
    println("Optimization Results:")
    println("Optimal x: ", sol.u)
    println("Optimal objective value: ", sol.objective)
    println("Optimized a: ", params.a)
    println("Optimized b: ", params.b)
    println("Solver status: ", sol.status)
end
