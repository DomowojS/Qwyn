#=
Optimisation module
=#

module Optimisation
include("Input_Processing.jl")
include("Initialisation_Module.jl")
include("SimpleComputation.jl")
using .SimpleComputation, .Input_Processing, .Initialisation_Module

export Test

function Test(x)
    
println("This Shit works")
println(x)

end



end