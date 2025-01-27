using NLopt, ForwardDiff, Revise


function autodiff(f::Function)
    function nlopt_fn(x::Vector, grad::Vector)
        if length(grad) > 0
            # Use ForwardDiff to compute the gradient. Replace with your
            # favorite Julia automatic differentiation package.
            ForwardDiff.gradient!(grad, f, x)
        end
        return f(x, P)
    end
end

function my_objective_fn(x::Vector, P)
    P.x=x
    return sqrt(P.x[2])
end

function my_constraint_fn(x::Vector, a, b)
    
     return (a * x[1] + b)^3 - x[2]
end

# These functions do not implement `grad`:
function Set_and_RunOptim()
    #Pre-Assign struct
    global P=params(0.0, [5.0, 1.0])

    opt = NLopt.Opt(:LD_MMA, 2)
    NLopt.lower_bounds!(opt, [-Inf, 0.0])
    NLopt.xtol_rel!(opt, 1e-4)
    # But we wrap them in autodiff before passing to NLopt:
    NLopt.min_objective!(opt, autodiff(x -> my_objective_fn(x, P)))
    NLopt.inequality_constraint!(opt, autodiff(x -> my_constraint_fn(x, 2, 0)), 1e-8)
    NLopt.inequality_constraint!(opt, autodiff(x -> my_constraint_fn(x, -1, 1)), 1e-8)
    min_f, min_x, ret = NLopt.optimize(opt, [1.234, 5.678])
    num_evals = NLopt.numevals(opt)
    println(
        """
        objective value       : $min_f
        solution              : $min_x
        solution status       : $ret
        # function evaluation : $num_evals
        """
    )
end

mutable struct params
    y::Float64
    x::Vector{Float64}
end