include("SmallestExample.jl")
using .SmallestExample

function Optim()
    P=params(0, [5.0, 1.0, 0.5, 0.5], 0)
    min_f, min_x, ret, num_evals = Set_and_RunOptim(P)

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
    n::Int64
end