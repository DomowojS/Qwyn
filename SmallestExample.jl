

module SmallestExample
include("FakeModule.jl")
using NLopt, .FakeModule, ForwardDiff
export Set_and_RunOptim, Optim

    function Optim(y, x)
        P=params(y, x)
        Set_and_RunOptim(P)
    end

    function my_constraint_fn(x::Vector, grad::Vector, a, b)
        return (a * x[1] + b)^3 - x[2]
    end

    function Set_and_RunOptim(P)
    trace=Any[]
    opt = NLopt.Opt(:LN_COBYLA, 2)
    NLopt.lower_bounds!(opt, [-Inf, 0.0])
    #NLopt.upper_bounds!(opt, [+Inf, 0.0])
    NLopt.xtol_rel!(opt, 1e-4)
    NLopt.min_objective!(opt, (x, grad) -> my_objective_fn(x, grad, P, trace))
    NLopt.inequality_constraint!(opt, (x, g) -> my_constraint_fn(x, g, 2, 0), 1e-8)
    NLopt.inequality_constraint!(opt, (x, g) -> my_constraint_fn(x, g, -1, 1), 1e-8)
    min_f, min_x, ret = NLopt.optimize(opt, P.x)
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

    function comp_of_obj(P)

        return sqrt(FakeFunct(P))

    end

    function my_objective_fn(x::Vector, grad::Vector, P, trace)
        P.x=x
        println("Evaluated")
        return comp_of_obj(P)
    end

    mutable struct params
        y::Float64
        x::Vector{Float64}
    end

end