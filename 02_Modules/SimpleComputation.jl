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
export Ishihara_WakeModel!, Superposition!, getTurbineInflow!, getNewThrustandPower!

# Compute single wake
function Ishihara_WakeModel!(WindFarm, CS)

    CS.k, CS.epsilon, CS.a, CS.b, CS.c, CS.d, CS.e, CS.f = ComputeEmpiricalVars(CS.Ct_vec, CS.TI_0_vec, 
                                                            CS.k, CS.epsilon, CS.a, CS.b, CS.c, CS.d, CS.e, CS.f); # Compute empirical values

    CS.Computation_Region_ID = (CS.XCoordinates .> 0.1e-10) .& (CS.YCoordinates .< 20 .* WindFarm.D) # Limit computation domain to reasonable scope

    # Representative wake width (sigma(x))
    CS.sigma    .= CS.Computation_Region_ID .* (CS.k .* CS.XCoordinates./WindFarm.D .+ CS.epsilon) .* WindFarm.D; # Compute wake width of all turbines
    
    # Velocity deficit
    CS.Delta_U  .=  CS.Computation_Region_ID .* ((1 ./ (CS.a .+ CS.b .* CS.XCoordinates./WindFarm.D .+ CS.c .* (1 .+ CS.XCoordinates./WindFarm.D).^-2).^2) .* exp.(-CS.r.^2 ./(2 .* CS.sigma.^2)) .* CS.u_0_vec);# Compute velocity deficit
    
    #Rotor-added turbulence
    #Include turbulence computation
    CS.k1       .=   (1 .- (CS.r./WindFarm.D .<= 0.5)) .+ (CS.r./WindFarm.D .<= 0.5) .* ((cos.(pi./2 .* (CS.r./WindFarm.D .- 0.5))).^2);
    CS.k2       .=   (CS.r./WindFarm.D .<= 0.5) .* ((cos.(pi./2 .* (CS.r./WindFarm.D .+ 0.5))).^2);
    CS.delta    .=   (CS.ZCoordinates .< WindFarm.H) .* (WindFarm.TI_a .* (sin.(pi .* (WindFarm.H .- (CS.ZCoordinates))./WindFarm.H)).^2);

    CS.Delta_TI .=   CS.Computation_Region_ID .* (((1 ./ (CS.d .+ CS.e .* CS.XCoordinates./WindFarm.D .+ CS.f .* (1 .+ CS.XCoordinates./WindFarm.D).^-2)) .* 
                            (CS.k1 .* exp.(-(CS.r .- 0.5.*WindFarm.D).^2 ./ (2 .* (CS.sigma).^2)) .+ CS.k2 .* exp.(-(CS.r .+ 0.5.*WindFarm.D).^2 ./(2 .* (CS.sigma).^2)))) .- CS.delta);# Compute rotor-added turbulence
end #Ishihara_Wakemodel
#

function ComputeEmpiricalVars(Ct, TI_0_vec, k, epsilon, a, b, c, d, e, f)
    k       .= 0.11 .* Ct.^1.07  .* TI_0_vec.^0.2 
    epsilon .= 0.23 .* Ct.^-0.25 .* TI_0_vec.^0.17
    a       .= 0.93 .* Ct.^-0.75 .* TI_0_vec.^0.17
    b       .= 0.42 .* Ct.^0.6   .* TI_0_vec.^0.2
    c       .= 0.15 .* Ct.^-0.25 .* TI_0_vec.^-0.7
    d       .= 2.3  .* Ct.^1.2   .* TI_0_vec.^0.1
    e       .= 1.0               .* TI_0_vec.^0.1
    f       .= 0.7  .* Ct.^-3.2  .* TI_0_vec.^-0.45
    return k, epsilon, a, b, c, d, e, f
end#ComputeEmpiricalVars

# Compute mixed wake properties
function Superposition!(WindFarm, CS)

    if WindFarm.Superpos == "Linear_Rotorbased"
        #Velocity deficit
        CS.U_Farm .= WindFarm.u_ambient_zprofile .- sum(CS.Delta_U, dims=3);
        CS.U_Farm[CS.U_Farm.<0].=0 #Filter of "negative" wind speeds at low heights

        #Rotor-added turbulence
        ### IMPLEMENT Height Profile for TI_a -> (WindFarm.TI_a.*WindFarm.u_ambient).^2 needs to be height related and ./WindFarm.u_ambient;, too!
        CS.TI_Farm .= sqrt.((WindFarm.TI_a.*WindFarm.u_ambient).^2 .+ sum((CS.Delta_TI.*CS.u_0_vec).^2, dims=3))./WindFarm.u_ambient;
    elseif WindFarm.Superpos == "Momentum_Conserving"
        CS.U_Farm .= CS.Delta_U;
    end
end#Superposition

#Evaluate new inflow data
function getTurbineInflow!(WindFarm, CS) 
    # Update inflow parameter
    CS.u_0_vec_old .= CS.u_0_vec #Store old inflow data
    CS.u_0_vec .= reshape(mean(mean(CS.U_Farm, dims=3), dims=2), (1,1,WindFarm.N))    #Compute mean inflow velocity for each turbine
    CS.TI_0_vec .= reshape(mean(mean(CS.TI_Farm, dims=3), dims=2), (1,1,WindFarm.N))  #Compute mean Turbulence intensity for each turbine


    # Compute termination criterion
    CS.zeta = findmax(abs.(CS.u_0_vec.-CS.u_0_vec_old))[1]
end#getTurbineInflow

# Compute new turbine properties
function getNewThrustandPower!(WindFarm, CS)
    if WindFarm.Turbine_Type=="VestasV80"
        CS.Ct_vec   .=  CS.Interp_Ct.(CS.u_0_vec);  #Ct of each turbine
        CS.P_vec    .=  CS.Interp_P.(CS.u_0_vec);   #P of each turbine 
    elseif WindFarm.Turbine_Type=="NREL_5MW"
        CS.Ct_vec   .=  CS.Interp_Ct.(CS.u_0_vec);  #Ct of each turbine
        CS.P_vec    .=  CS.Interp_P.(CS.u_0_vec);   #P of each turbine
    elseif WindFarm.Turbine_Type=="DTU_10MW"
        CS.Ct_vec   .=  CS.Interp_Ct.(CS.u_0_vec);  #Ct of each turbine
        CS.P_vec    .=  CS.Interp_P.(CS.u_0_vec);   #P of each turbine
    elseif WindFarm.Turbine_Type=="IEA_15MW"
        CS.Ct_vec   .=  CS.Interp_Ct.(CS.u_0_vec);  #Ct of each turbine
        CS.P_vec    .=  CS.Interp_P.(CS.u_0_vec);   #P of each turbine
    end
end

# Compute power output



# Functions for single wake computation

end