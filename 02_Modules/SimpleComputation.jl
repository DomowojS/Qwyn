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
using Statistics
export Ishihara_WakeModel!, Superposition!, getTurbineInflow!

# Compute single wake
function Ishihara_WakeModel!(WindFarm, CS)
    println("computing single wake..")
    println("..computing empirical values..")
    CS.k, CS.epsilon, CS.a, CS.b, CS.c, CS.d, CS.e, CS.f = ComputeEmpiricalVars(CS.Ct_vec, WindFarm.TI_a, 
                                                            CS.k, CS.epsilon, CS.a, CS.b, CS.c, CS.d, CS.e, CS.f); # Compute empirical values
    println("..computing r..")
    CS.r        .= sqrt.((CS.YCoordinates.*WindFarm.D).^2 .+ (CS.Z_Levels.-WindFarm.H).^2) # Compute vector in radial & height direction for computation

    CS.Computation_Region_ID = (CS.XCoordinates .> 0.1e-10) .& (CS.YCoordinates .< 20) # Limit computation domain to reasonable scope

    # Representative wake width (sigma(x))
    println("..computing sigma..")
    CS.sigma    .= CS.Computation_Region_ID .* (CS.k .* CS.XCoordinates .+ CS.epsilon) .* WindFarm.D; # Compute wake width of all turbines
    
    # Velocity deficit
    println("..computing velocity deficit..")
    CS.Delta_U  .=  CS.Computation_Region_ID .* ((1 ./ (CS.a .+ CS.b .* CS.XCoordinates .+ CS.c .* (1 .+ CS.XCoordinates).^-2).^2) .* exp.(-CS.r.^2 ./(2 .* CS.sigma.^2)) .* CS.u_0_vec);# Compute velocity deficit
    
    #Rotor-added turbulence
    println("..computing turbulence empirical values..")
    #Include turbulence computation
    CS.k1       .=   (1 .- (CS.r .<= 0.5 * WindFarm.D)) .+ (CS.r .<= 0.5 * WindFarm.D) .* ((cos.(pi./2 .* (CS.r./WindFarm.D .- 0.5))).^2);
    CS.k2       .=   (CS.r .<= 0.5 * WindFarm.D) .* ((cos.(pi./2 .* (CS.r./WindFarm.D .+ 0.5))).^2);
    CS.delta    .=   (CS.Z_Levels .< WindFarm.H) .* (WindFarm.TI_a .* (sin.(pi .* (WindFarm.H .- (CS.Z_Levels))./WindFarm.H)).^2);

    println("..computing rotor-added turbulence..")
    CS.Delta_TI .=   CS.Computation_Region_ID .* (((1 ./ (CS.d .+ CS.e .* CS.XCoordinates .+ CS.f .* (1 .+ CS.XCoordinates).^-2)) .* 
                            (CS.k1 .* exp.(-(CS.r .- 0.5.*WindFarm.D).^2 ./ (2 .* (CS.sigma).^2)) .+ CS.k2 .* exp.(-(CS.r .+ 0.5.*WindFarm.D).^2 ./(2 .* (CS.sigma).^2)))) .- CS.delta);# Compute rotor-added turbulence
    

    println("..single wake computation finished!")
end #Ishihara_Wakemodel
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
end#ComputeEmpiricalVars

# Compute mixed wake properties
function Superposition!(WindFarm, CS)

    if WindFarm.Superpos == "Quadratic_Rotorbased"
        println("computing mixed wake region..")
        #Velocity deficit
        println("..computing velocity deficit..")
        CS.U_Farm .= WindFarm.u_ambient_zprofile .- sqrt.(sum(CS.Delta_U.^2, dims=4)); 
        CS.U_Farm[CS.U_Farm.<0].=0 #Filter of "negative" wind speeds in low heights

        #Rotor-added turbulence
        ### IMPLEMENT Height Profile for TI_a -> (WindFarm.TI_a.*WindFarm.u_ambient).^2 needs to be height related and ./WindFarm.u_ambient;, too!
        println("..computing rotor-added turbulence..")
        CS.TI_Farm .= sqrt.((WindFarm.TI_a.*WindFarm.u_ambient).^2 .+ sum((CS.Delta_TI.*CS.u_0_vec).^2, dims=4))./WindFarm.u_ambient;

    elseif WindFarm.Superpos =="Momentum_Conserving"
        CS.U_Farm .= CS.Delta_U;
    end

 println("..mixed wake computation finished!")
end#Superposition

#Evaluate new inflow data
function getTurbineInflow!(WindFarm, CS) 
    println("Reevaluating inflow parameter..")
    # Update inflow parameter
    CS.u_0_vec_old .= CS.u_0_vec #Store old inflow data
    CS.u_0_vec .= reshape(mean(mean(CS.U_Farm, dims=3), dims=2), (1,1,1,WindFarm.N))    #Compute mean inflow velocity for each turbine
    CS.TI_0_vec .= reshape(mean(mean(CS.TI_Farm, dims=3), dims=2), (1,1,1,WindFarm.N))  #Compute mean Turbulence intensity for each turbine

    # Compute termination criterion
    CS.zeta = findmax(abs.(CS.u_0_vec.-CS.u_0_vec_old))[1]
    println("finished reevaluation!")
    
# Convert the struct to a dictionary
struct_dict = Dict{String, Any}(string.(propertynames(CS)) .=> getfield.(Ref(CS), propertynames(CS)))
struct_dict2 = Dict{String, Any}(string.(propertynames(WindFarm)) .=> getfield.(Ref(WindFarm), propertynames(WindFarm)))

# Specify the filename for the .mat file
filename = "99_PlotWMATLAB/WindFarmCS_Test.mat"
filename2 = "99_PlotWMATLAB/WindFarmWF.mat"
# Save the struct to the .mat file
matwrite(filename, struct_dict)
matwrite(filename2, struct_dict2)
end#getTurbineInflow

# Compute new turbine properties

# Compute power output



# Functions for single wake computation

end