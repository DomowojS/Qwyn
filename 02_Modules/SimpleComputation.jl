#= 
Script to Process computation &/ optimisation
=#

module SimpleComputation
using Statistics
export Ishihara_WakeModel!, Superposition!, getTurbineInflow!, getNewThrustandPower!, getTotalPower!

function Ishihara_WakeModel!(WindFarm, CS)
# Compute single wake according to the Ishihara-Qian model (2018)
    CS.k, CS.epsilon, CS.a, CS.b, CS.c, CS.d, CS.e, CS.f = ComputeEmpiricalVars(CS.Ct_vec, CS.TI_0_vec, 
                                                            CS.k, CS.epsilon, CS.a, CS.b, CS.c, CS.d, CS.e, CS.f); # Compute empirical values

    CS.Computation_Region_ID = (CS.XCoordinates .> 0.1e-10) .& (CS.YCoordinates .< 20 .* WindFarm.D) # Limit computation domain to reasonable scope

    # Representative wake width (sigma(x))
    CS.sigma .= CS.Computation_Region_ID .* (CS.k .* CS.XCoordinates./WindFarm.D .+ CS.epsilon) .* WindFarm.D; # Compute wake width of all turbines
    
    # Compute correction terms & convection velocity if needed
    if WindFarm.Meandering == true 	|| WindFarm.Superpos == "Momentum_Conserving"

        CS.u_c_vec = comp_ConvectionVelocity(CS.tmp, CS.u_c_vec, CS.Ct_vec, WindFarm.D, CS.sigma, CS.u_0_vec)

        # Only for meandering model
        if WindFarm.Meandering == true
        # run meandering correction
        CS.psi, CS.Lambda, CS.sigma_m = Meandering_Correction(CS.sigma_m, CS.psi, CS.Lambda, CS.TI_0_vec, CS.u_0_vec, CS.ZCoordinates, CS.XCoordinates, CS.u_c_vec)
        # Corrected velocity deficit
        CS.Delta_U  .=  CS.Computation_Region_ID .* ((1 ./ (CS.a .+ CS.b .* CS.XCoordinates./WindFarm.D .+ CS.c .* (1 .+ CS.XCoordinates./WindFarm.D).^-2).^2) .* (1 .+ (CS.sigma_m./CS.sigma).^2).^-0.5 .* exp.(-CS.r.^2 ./(2 .* (CS.sigma.^2 .+ CS.sigma_m.^2))) .* CS.u_0_vec);# Compute velocity deficit
    	else
        # Velocity deficit without correction
        CS.Delta_U  .=  CS.Computation_Region_ID .* ((1 ./ (CS.a .+ CS.b .* CS.XCoordinates./WindFarm.D .+ CS.c .* (1 .+ CS.XCoordinates./WindFarm.D).^-2).^2) .* exp.(-CS.r.^2 ./(2 .* CS.sigma.^2)) .* CS.u_0_vec);# Compute velocity deficit
        end
    else

        # Velocity deficit without corrections
        CS.Delta_U  .=  CS.Computation_Region_ID .* ((1 ./ (CS.a .+ CS.b .* CS.XCoordinates./WindFarm.D .+ CS.c .* (1 .+ CS.XCoordinates./WindFarm.D).^-2).^2) .* exp.(-CS.r.^2 ./(2 .* CS.sigma.^2)) .* CS.u_0_vec);# Compute velocity deficit
    end

    #Rotor-added turbulence
    #Include turbulence computation
    CS.k1       .=   (1 .- (CS.r./WindFarm.D .<= 0.5)) .+ (CS.r./WindFarm.D .<= 0.5) .* ((cos.(pi./2 .* (CS.r./WindFarm.D .- 0.5))).^2);
    CS.k2       .=   (CS.r./WindFarm.D .<= 0.5) .* ((cos.(pi./2 .* (CS.r./WindFarm.D .+ 0.5))).^2);
    CS.delta    .=   (CS.ZCoordinates .< WindFarm.H) .* (WindFarm.TI_a .* (sin.(pi .* (WindFarm.H .- (CS.ZCoordinates))./WindFarm.H)).^2);

    CS.Delta_TI .=   CS.Computation_Region_ID .* (((1 ./ (CS.d .+ CS.e .* CS.XCoordinates./WindFarm.D .+ CS.f .* (1 .+ CS.XCoordinates./WindFarm.D).^-2)) .* 
                            (CS.k1 .* exp.(-(CS.r .- 0.5.*WindFarm.D).^2 ./ (2 .* (CS.sigma).^2)) .+ CS.k2 .* exp.(-(CS.r .+ 0.5.*WindFarm.D).^2 ./(2 .* (CS.sigma).^2)))) .- CS.delta);# Compute rotor-added turbulence
end #Ishihara_Wakemodel
#

function ComputeEmpiricalVars(Ct::Array{Float64}, TI_0_vec::Array{Float64}, k::Array{Float64}, epsilon::Array{Float64}, a::Array{Float64}, b::Array{Float64}, c::Array{Float64}, d::Array{Float64}, e::Array{Float64}, f::Array{Float64})
# Compute empirical parameters for the Ishihara WakeModel
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

function Meandering_Correction(sigma_m::Array{Float64}, psi::Array{Float64}, Lambda::Array{Float64}, TI_0_vec::Array{Float64}, u_0_vec::Array{Float64}, ZCoordinates::Array{Float64}, XCoordinates::Array{Float64}, u_c_vec::Array{Float64})
# Compute the meandering correction according to the Braunbehrens & Segalini model (2019) 
    #Compute fluctuation intensity
    psi     .= 0.7.*TI_0_vec.*u_0_vec    
    #Compute integral length scale of representative eddy   
    Lambda  .= (0.4 .* ZCoordinates) ./ psi
    #Compute corrected wake width 
    sigma_m .= sqrt.((2 .* psi .* Lambda.^2) .* ((XCoordinates./(u_c_vec.*Lambda)) .+ exp.(-XCoordinates./(u_c_vec.*Lambda)) .- 1))  
    return psi, Lambda, sigma_m          
end#Meandering_Correction

function comp_ConvectionVelocity(tmp, u_c_vec, Ct, D, sigma, u_0_vec)
# Compute local convection velocity of each turbine
    # Term within the sqrt. expression
    tmp = (1 .- ((Ct.*D^2)./(8 .* sigma.^2)))
    tmp[tmp.==-Inf] .= 0 #Correction (-Inf appears for those points which should not be computed, since they are dowmstream)
    u_c_vec .= (0.5 .+ 0.5 .* sqrt.(tmp)) .* u_0_vec
    return u_c_vec
end#comp_ConvectionVelocity

function Superposition!(WindFarm, CS)
# Compute mixed wake properties
    if WindFarm.Superpos == "Linear_Rotorbased"
        #Velocity deficit
        CS.U_Farm .= WindFarm.u_ambient_zprofile .- sum(CS.Delta_U, dims=3);
        CS.U_Farm[CS.U_Farm.<0].=0 #Filter of "negative" wind speeds at low heights

        #Rotor-added turbulence
        ### IMPLEMENT Height Profile for TI_a -> (WindFarm.TI_a.*WindFarm.u_ambient).^2 needs to be height related and ./WindFarm.u_ambient;, too!
        CS.TI_Farm .= sqrt.((WindFarm.TI_a.*WindFarm.u_ambient).^2 .+ sum((CS.Delta_TI.*CS.u_0_vec).^2, dims=3))./WindFarm.u_ambient;
    elseif WindFarm.Superpos == "Momentum_Conserving"
        #Compute global wake convection velocity
        if CS.i == 1 #For first iteration U_c_Farm gets initial conditions
            CS.U_c_Farm=CS.u_c_vec;

        else

            CS.U_Farm .= CS.Delta_U;
    
        end
        
        
    end
end#Superposition


function getTurbineInflow!(WindFarm, CS) 
# Evaluate new inflow data
    CS.u_0_vec_old .= CS.u_0_vec #Store old inflow data
    CS.u_0_vec .= reshape(mean(mean(CS.U_Farm, dims=3), dims=2), (1,1,WindFarm.N))    #Compute mean inflow velocity for each turbine
    CS.TI_0_vec .= reshape(mean(mean(CS.TI_Farm, dims=3), dims=2), (1,1,WindFarm.N))  #Compute mean Turbulence intensity for each turbine


    # Compute termination criterion
    CS.zeta = findmax(abs.(CS.u_0_vec.-CS.u_0_vec_old))[1]
end#getTurbineInflow


function getNewThrustandPower!(WindFarm, CS)
# Compute new turbine properties
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

function getTotalPower!(CS)
# Compute total pwoer output of the wind farm
    CS.TotalPower = sum(CS.P_vec)
end#getTotalPower

# Compute power output



# Functions for single wake computation

end