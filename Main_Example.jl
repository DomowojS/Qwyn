include("SmallestExample.jl")
using .SmallestExample

function Optim()
    P=params(0, [5.0, 1.0])
    Set_and_RunOptim(P)
end

mutable struct params
    y::Float64
    x::Vector{Float64}
end