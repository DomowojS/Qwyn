using JuMP
using Ipopt
using LinearAlgebra  # For distance calculations

# Mutable struct containing x and y vectors
mutable struct Parameters
    x::Vector{Float64}
    y::Vector{Float64}
    a::Float64
    b::Float64
    c::Float64
end

# Your objective function
function AEP4Optim(params)
    # Replace this with your actual function
    # Now using x and y from params
    return sum(params.x.^2) + params.a * sum(params.y.^2) + params.b * sum(params.x.*params.y) + params.c
end

# Calculate distance between two points
function distance(x1, y1, x2, y2)
    return sqrt((x1 - x2)^2 + (y1 - y2)^2)
end

# Constraint function for vector x
function x_constraint(y)
    return y.^2  # Element-wise operation for vectors
end

# Your objective function that takes vectors and parameters
function AEP4Optim(x, y, params)
    # Replace this with your actual function
    # This is just a dummy example
    return sum(x.^2) + params.a * sum(y.^2) + params.b * sum(x.*y) + params.c
end

# Main optimization
function optimize_layout()
    n = length(params.x)  # Get length from existing vectors
    min_distance = 3.7    # Minimum distance between machines
    
    # Create model
    model = Model(Ipopt.Optimizer)

    # Define variables
    @variable(model, x[1:n])
    @variable(model, y[1:n])

    # Set start values
    for i in 1:n
        set_start_value(x[i], params.x[i])
        set_start_value(y[i], params.y[i])
    end

    # Set any bounds if needed
    for i in 1:n
        set_lower_bound(y[i], -15)
        set_upper_bound(y[i], 15)
    end

    # Add minimum distance constraints between all pairs of machines
    for i in 1:n
        for j in (i+1):n
            @NLconstraint(model, 
                sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2) >= min_distance
            )
        end
    end

    # Add your x constraints if still needed
    for i in 1:n
        @constraint(model, x[i] == x_constraint(y[i]))
    end

    # Create modified parameters struct for objective function
    # that uses the optimization variables
    @objective(model, Max, AEP4Optim(x, y, params))

    # Solve
    optimize!(model)

    # Get results
    x_optimal = value.(x)
    y_optimal = value.(y)
    optimal_value = objective_value(model)

    return x_optimal, y_optimal, optimal_value
end

# Example usage:
# Assume params is already defined with initial values
params = Parameters(
    rand(5),  # Initial x values
    rand(5),  # Initial y values
    2.0,      # a
    1.5,      # b
    0.5       # c
)

x_optimal, y_optimal, final_value = optimize_layout()

# Print results
println("Optimal value: ", final_value)
println("\nInitial positions:")
for i in 1:length(params.x)
    println("Machine $i: (", params.x[i], ", ", params.y[i], ")")
end
println("\nOptimal positions:")
for i in 1:length(x_optimal)
    println("Machine $i: (", x_optimal[i], ", ", y_optimal[i], ")")
end

# Verify minimum distances are maintained
println("\nChecking minimum distances:")
for i in 1:length(x_optimal)
    for j in (i+1):length(x_optimal)
        d = distance(x_optimal[i], y_optimal[i], x_optimal[j], y_optimal[j])
        println("Distance between machines $i and $j: $d")
    end
end