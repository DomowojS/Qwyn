module SmallestExample
include("FakeModule.jl")
using NLopt, .FakeModule, ForwardDiff
export Set_and_RunOptim, Optim

    function Set_and_RunOptim(P)
        trace = Any[]
        N = length(P.x) ÷ 2  # Number of machines
        opt = NLopt.Opt(:LN_COBYLA, length(P.x))
        
        # Set bounds for all coordinates
        NLopt.lower_bounds!(opt, fill(-Inf, length(P.x)))
        NLopt.upper_bounds!(opt, fill(+Inf, length(P.x)))
        NLopt.xtol_rel!(opt, 1e-4)
        
        NLopt.min_objective!(opt, (x, grad) -> my_objective_fn(x, grad, P, trace))
        
        # Original scalar constraints
        NLopt.inequality_constraint!(opt, (x, g) -> my_constraint_fn(x, g, 2, 0), 1e-8)
        NLopt.inequality_constraint!(opt, (x, g) -> my_constraint_fn(x, g, -1, 1), 1e-8)
        
        # Add pairwise distance constraints
        for i in 1:N
            for j in (i+1):N
                NLopt.inequality_constraint!(opt, 
                    (x, g) -> distance_constraint(x, g, i, j, 1.0), 
                    1e-8)
            end
        end
        
        min_f, min_x, ret = NLopt.optimize(opt, P.x)
        num_evals = NLopt.numevals(opt)

        return min_f, min_x, ret, num_evals
    end

    function comp_of_obj(P, x)
        return sqrt(FakeFunct(P, x))
    end

    function my_objective_fn(x::Vector, grad::Vector, P, trace)
        println("Evaluated")
        return comp_of_obj(P, x)
    end

    function my_constraint_fn(x::Vector, grad::Vector, a, b)
        return (a * x[1] + b)^3 - x[2]
    end

    # Helper function to compute distance between two machines
    function get_machine_distance(x::Vector, i::Int, j::Int)
        N = length(x) ÷ 2
        # Get coordinates for machine i
        xi, yi = x[i], x[i + N]
        # Get coordinates for machine j
        xj, yj = x[j], x[j + N]
        # Compute Euclidean distance
        return sqrt((xi - xj)^2 + (yi - yj)^2)
    end

    # Distance constraint function
    # Returns positive value if constraint is violated (distance > max_dist)
    # Returns negative value if constraint is satisfied (distance < max_dist)
    function distance_constraint(x::Vector, grad::Vector, i::Int, j::Int, max_dist::Float64)
        dist = get_machine_distance(x, i, j)
        return max_dist - dist
    end
end