#= 
Script to Process computation &/ optimisation
=#

module SimpleComputation
using Pkg
using FileIO
using Parameters
using InteractiveUtils
using DataStructures
using LatinHypercubeSampling
using Revise
#using PlotlyJS
using MAT
using LinearAlgebra
export Ishihara_WakeModel

function Ishihara_WakeModel(WindFarm, CS)
    println("computing single wake..")
    println("..computing empirical values..")
    CS.k, CS.epsilon, CS.a, CS.b, CS.c, CS.d, CS.e, CS.f = ComputeEmpiricalVars(CS.Ct_vec, WindFarm.TI_a, 
                                                            CS.k, CS.epsilon, CS.a, CS.b, CS.c, CS.d, CS.e, CS.f); # Compute empirical values
    println("..computing r..")
    CS.r = sqrt.((CS.YCoordinates.*WindFarm.D).^2 .+ ((CS.Z_Levels.*WindFarm.D).-WindFarm.H).^2) # Compute vector in radial & height direction for computation

    # Representative wake width (sigma(x))
    println("..computing sigma..")
    CS.sigma    =   ifelse.((CS.XCoordinates .> 0.1e-10) .& (CS.YCoordinates .< 20), (CS.k .* CS.XCoordinates .+ CS.epsilon) .* WindFarm.D, 0); # Compute wake width of all turbines
    # Velocity deficit
    println("..computing velocity deficit..")
    CS.Delta_U   =  ifelse.((CS.XCoordinates .> 0.1e-10) .& (CS.YCoordinates .< 20), (1 ./ (CS.a .+ CS.b .* CS.XCoordinates .+ CS.c .* (1 .+ CS.XCoordinates).^-2).^2) .* exp.(-CS.r.^2 ./(2 .* CS.sigma.^2)) .* CS.u_0_vec, 0);# Compute velocity deficit
    
    #Rotor-added turbulence
    println("..computing turbulence empirical values..")
    #Include turbulence computation
    CS.k1       =   ifelse.(CS.r .<= 0.5 * WindFarm.D, (cos.(pi./2 .* (CS.r./WindFarm.D .- 0.5))).^2, 1);
    CS.k2       =   ifelse.(CS.r .<= 0.5 * WindFarm.D, (cos.(pi./2 .* (CS.r./WindFarm.D .+ 0.5))).^2, 0);
    CS.delta    =   ifelse.(CS.Z_Levels.*WindFarm.D .< WindFarm.H, WindFarm.TI_a .* (sin.(pi .* (WindFarm.H .- (CS.Z_Levels.*WindFarm.D))./WindFarm.H)).^2, 0);

    println("..computing rotor-added turbulence..")
    CS.Delta_TI =   ifelse.((CS.XCoordinates .> 0.1e-10) .& (CS.YCoordinates .< 20), ((1 ./ (CS.d .+ CS.e .* CS.XCoordinates .+ CS.f .* (1 .+ CS.XCoordinates).^-2)) .* 
                            (CS.k1 .* exp.(-(CS.r .- 0.5.*WindFarm.D).^2 ./(2 .* (CS.sigma).^2)) .+ CS.k2 .* exp.(-(CS.r .+ 0.5.*WindFarm.D).^2 ./(2 .* (CS.sigma).^2)))) .- CS.delta,
                            0);# Compute rotor-added turbulence
    

    println("..single wake computation finished!")
                            
                            #
# Convert the struct to a dictionary
struct_dict = Dict{String, Any}(string.(propertynames(CS)) .=> getfield.(Ref(CS), propertynames(CS)))
struct_dict2 = Dict{String, Any}(string.(propertynames(WindFarm)) .=> getfield.(Ref(WindFarm), propertynames(WindFarm)))

# Specify the filename for the .mat file
filename = "99_PlotWMATLAB/WindFarmCS_CheckNewCoordinates.mat"
filename2 = "99_PlotWMATLAB/WindFarmWF_CheckNewCoordinates.mat"
# Save the struct to the .mat file
matwrite(filename, struct_dict)
matwrite(filename2, struct_dict2)

    return WindFarm, CS
end
#

function ComputeEmpiricalVars(Ct, TI_a, k, epsilon, a, b, c, d, e, f)
    k       .= 0.11 .* Ct.^1.07  .* TI_a^0.2 
    epsilon .= 0.23 .* Ct.^-0.25 .* TI_a^0.17
    a       .= 0.93 .* Ct.^-0.75 .* TI_a^0.17
    b       .= 0.42 .* Ct.^0.6   .* TI_a^0.2
    c       .= 0.15 .* Ct.^-0.25 .* TI_a^-0.7
    d       .= 2.3  .* Ct.^1.2   .* TI_a^0.1
    e       .= 1.0               .* TI_a^0.1
    f       .= 0.7  .* Ct.^-3.2  .* TI_a^-0.45
    return k, epsilon, a, b, c, d, e, f
end
# Compute single wake

# Compute mixed wake 

# Evaluate new inflow data

# Compute new turbine properties

# Compute power output



# Functions for single wake computation

end