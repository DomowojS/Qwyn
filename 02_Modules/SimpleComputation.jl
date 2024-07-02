#= 
Script to Process computation &/ optimisation
=#

module SimpleComputation
using Statistics
export FindComputationRegion!, Single_Wake_Computation!, Superposition!, getTurbineInflow!, getNewThrustandPower!, getTotalPower!, computeTerminationCriterion!

function FindComputationRegion!(WindFarm, CS)
# Find computation region (to prevent unneccasary computation)
    CS.Computation_Region_ID        = (CS.XCoordinates .> 0.1e-10) .& (abs.(CS.YCoordinates) .< 20 .* WindFarm.D)  # Limit computation domain to reasonable scope
    CS.Computation_Region_ID_for_Uc  = (CS.XCoordinates .> 0.1e-10) .& (abs.(CS.Y_for_Uc) .< 20 .* WindFarm.D)
    CS.ID_Turbines              = vec((CS.Ct_vec .> 0) .& (CS.zetaID .!= 0))  #Identify turbines which should not be included in computation
    CS.ID_OutOperConst          = vec(CS.Ct_vec .== 0)                     #Identify only turbines which are outside of operation constraints (below cut in/ above cut out windspeed)    
end

function Single_Wake_Computation!(WindFarm, CS)
## For Ishihara - single wake model:
    CS.k, CS.epsilon, CS.a, CS.b, CS.c, CS.d, CS.e, CS.f = ComputeEmpiricalVars(CS.Ct_vec, CS.TI_0_vec, 
    CS.k, CS.epsilon, CS.a, CS.b, CS.c, CS.d, CS.e, CS.f, CS.ID_Turbines); # Compute empirical values
    
    Ishihara_WakeModel!(WindFarm, CS, CS.Delta_U, CS.Delta_TI, CS.XCoordinates, CS.ZCoordinates, CS.r, CS.Computation_Region_ID, CS.sigma, CS.sigma_m, CS.Lambda, CS.k1, CS.k2, CS.delta, false)   
    if WindFarm.Superpos == "Momentum_Conserving"
        Ishihara_WakeModel!(WindFarm, CS, CS.Delta_U_for_Uc, CS. Delta_TI_for_Uc, CS.XCoordinates, CS.Z_for_Uc, CS.r_for_Uc, CS.Computation_Region_ID_for_Uc, CS.sigma_for_Uc, CS.sigma_m_for_Uc, CS.Lambda_for_Uc, CS.k1_for_Uc, CS.k2_for_Uc, CS.delta_for_Uc, true) 
    end

end#Single_Wake_Computation

function Superposition!(WindFarm, CS)
# Compute mixed wake properties
    #Velocity deficit
    if WindFarm.Superpos == "Linear_Rotorbased"
        #Compute linear rotorbased sum
        CS.U_Farm .= WindFarm.u_ambient_zprofile .- sum(CS.Delta_U, dims=3);

    elseif WindFarm.Superpos == "Momentum_Conserving"
        
        #For first iteration U_c_Farm gets initial conditions
        CS.U_c_Farm = maximum(CS.u_c_vec, dims=3);
        CS.U_c_Farm[CS.U_c_Farm .== 0] .= 1;         #Correction, to prevent NaN in next computation step.
        #Compute weighting factor
        CS.weighting_Factor_for_Uc .= (CS.u_c_vec./CS.U_c_Farm)
        CS.weighting_Factor_for_Uc[CS.Delta_U_for_Uc .< (0.1 .* CS.u_ambient_for_Uc)] .= 1 #correction. For u_i < 0.1 of ambient wind speed -> no weighting is considered.
        #Compute weighted sum
        CS.Mixed_wake_for_Uc .= CS.u_ambient_for_Uc .- sum((CS.weighting_Factor_for_Uc .* CS.Delta_U_for_Uc), dims=3);
        
        #Compute global wake convection velocity
        i=0
        while any(abs.((CS.U_c_Farm_old .- CS.U_c_Farm) ./ CS.U_c_Farm) .>= 0.001) == true #any(abs.((CS.U_c_Farm_old .- CS.U_c_Farm) ./ CS.U_c_Farm) .>= 0.001) == true
            i=i+1
            #Safe old convection velocity for termination criterion computation
            CS.U_c_Farm_old .= CS.U_c_Farm
            #For iteration >2 compute U_c_Farm from last iteration's result
            CS.U_c_Farm .= sum((CS.Mixed_wake_for_Uc .* (CS.u_ambient_for_Uc .- CS.Mixed_wake_for_Uc)), dims=2) ./ sum((CS.u_ambient_for_Uc .- CS.Mixed_wake_for_Uc), dims=2);
            #CS.U_c_Farm[isnan.(CS.U_c_Farm)] .= WindFarm.u_ambient 
            #Compute weighting factor
            CS.weighting_Factor_for_Uc .= (CS.u_c_vec./CS.U_c_Farm)
            CS.weighting_Factor_for_Uc[CS.Delta_U_for_Uc .< (0.1 .* CS.u_ambient_for_Uc)] .= 1  #correction. For u_i < 0.1 of ambient wind speed -> no weighting is considered.
            #Compute weighted sum
            CS.Mixed_wake_for_Uc .= CS.u_ambient_for_Uc .- sum((CS.weighting_Factor_for_Uc .* CS.Delta_U_for_Uc), dims=3);

            println("Superpos-iteration: ", i)
        end

        # Now Compute the wake for the rotor points with converged farm convection velocity U_c_Farm
        CS.weighting_Factor .= (CS.u_c_vec./CS.U_c_Farm)
        CS.weighting_Factor[CS.Delta_U .< (0.1 .* WindFarm.u_ambient_zprofile)] .= 1 #correction. For u_i < 0.1 of ambient wind speed -> no weighting is considered.
        CS.U_Farm .= WindFarm.u_ambient_zprofile .- sum((CS.weighting_Factor .* CS.Delta_U), dims=3); 
              
    else 
        error("Wrong choice of superposition method. Check 'Superpos' input. Possible entries: 'Linear_Rotorbased' and 'Momentum_Conserving'.")   
    end
    CS.U_Farm[CS.U_Farm.<0].=0 #Filter of "negative" wind speeds at low heights

    #Rotor-added turbulence
    ### IMPLEMENT Height Profile for TI_a -> (WindFarm.TI_a.*WindFarm.u_ambient).^2 needs to be height related and ./WindFarm.u_ambient;, too!
    CS.TI_Farm .= sqrt.((WindFarm.TI_a.*WindFarm.u_ambient).^2 .+ sum((CS.Delta_TI.*CS.u_0_vec).^2, dims=3))./WindFarm.u_ambient;
end#Superposition

function getTurbineInflow!(WindFarm, CS) 
# Evaluate new inflow data
    CS.u_0_vec_old .= CS.u_0_vec #Store old inflow data
    CS.u_0_vec .= reshape(mean(mean(CS.U_Farm, dims=3), dims=2), (1,1,WindFarm.N))    #Compute mean inflow velocity for each turbine
    CS.TI_0_vec .= reshape(mean(mean(CS.TI_Farm, dims=3), dims=2), (1,1,WindFarm.N))  #Compute mean Turbulence intensity for each turbine
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
end#getTotalPower#

function computeTerminationCriterion!(CS)
# Compute termination criterion zeta
    CS.zetaID   .= abs.(CS.u_0_vec.-CS.u_0_vec_old)
    CS.zeta     = findmax(CS.zetaID)[1]
    
end#computeTerminationCriterion


##### Subfunctions for single wake computation #######

function Ishihara_WakeModel!(WindFarm, CS, Delta_U, Delta_TI, X, Z, R, ID, sigma, sigma_m, Lambda, k1, k2, delta, Call_for_Uc)
# Compute single wake according to the Ishihara-Qian model (2018)

        # Representative wake width (sigma(x))
        sigma .= ID .* (CS.k .* X./WindFarm.D .+ CS.epsilon) .* WindFarm.D; # Compute wake width of all turbines

        # Compute correction terms & convection velocity if needed
        if WindFarm.Meandering == true || WindFarm.Superpos == "Momentum_Conserving"
            
            if Call_for_Uc ==false # Only computed during first call (for rotor region computation)
            CS.u_c_vec = comp_ConvectionVelocity(CS.tmp, CS.u_c_vec, CS.Ct_vec, WindFarm.D, sigma, CS.u_0_vec)
            end

            # Only for meandering model
            if WindFarm.Meandering == true
            # run meandering correction
            CS.psi, Lambda, sigma_m = Meandering_Correction(sigma_m, CS.psi, Lambda, CS.TI_0_vec, CS.u_0_vec, Z, X, CS.u_c_vec, Call_for_Uc)
            # Corrected velocity deficit
            Delta_U  .=  ID .* ((1 ./ (CS.a .+ CS.b .* X./WindFarm.D .+ CS.c .* (1 .+ X./WindFarm.D).^-2).^2) .* (1 .+ (sigma_m./sigma).^2).^-0.5 .* exp.(-R.^2 ./(2 .* (sigma.^2 .+ sigma_m.^2))) .* CS.u_0_vec);# Compute velocity deficit
            else
            # Velocity deficit without correction
            Delta_U  .=  ID .* ((1 ./ (CS.a .+ CS.b .* X./WindFarm.D .+ CS.c .* (1 .+ X./WindFarm.D).^-2).^2) .* exp.(-R.^2 ./(2 .* sigma.^2)) .* CS.u_0_vec);# Compute velocity deficit
            end
        else
    
            # Velocity deficit without correction / computation of convection velocity
            Delta_U  .=  ID .* ((1 ./ (CS.a .+ CS.b .* X./WindFarm.D .+ CS.c .* (1 .+ X./WindFarm.D).^-2).^2) .* exp.(-R.^2 ./(2 .* sigma.^2)) .* CS.u_0_vec);# Compute velocity deficit
        end
    
        #Rotor-added turbulence
        #Include turbulence computation
        k1       .=   (1 .- (R./WindFarm.D .<= 0.5)) .+ (R./WindFarm.D .<= 0.5) .* ((cos.(pi./2 .* (R./WindFarm.D .- 0.5))).^2);
        k2       .=   (R./WindFarm.D .<= 0.5) .* ((cos.(pi./2 .* (R./WindFarm.D .+ 0.5))).^2);
        delta    .=   (Z .< WindFarm.H) .* (WindFarm.TI_a .* (sin.(pi .* (WindFarm.H .- (Z))./WindFarm.H)).^2);
    
         Delta_TI .=   ID .* (((1 ./ (CS.d .+ CS.e .* X./WindFarm.D .+ CS.f .* (1 .+ X./WindFarm.D).^-2)) .* 
                                (k1 .* exp.(-(R .- 0.5.*WindFarm.D).^2 ./ (2 .* (sigma).^2)) .+ k2 .* exp.(-(R .+ 0.5.*WindFarm.D).^2 ./(2 .* (sigma).^2)))) .- delta);# Compute rotor-added turbulence
end#Ishihara_Wakemodel


function ComputeEmpiricalVars(Ct::Array{Float64}, TI_0_vec::Array{Float64}, k::Array{Float64}, epsilon::Array{Float64}, a::Array{Float64}, b::Array{Float64}, c::Array{Float64}, d::Array{Float64}, e::Array{Float64}, f::Array{Float64}, ID_Turbines::BitVector)
# Compute empirical parameters for the Ishihara WakeModel
    k[:,:,ID_Turbines]       .= 0.11 .* Ct[:,:,ID_Turbines].^1.07  .* TI_0_vec[:,:,ID_Turbines].^0.2 
    epsilon[:,:,ID_Turbines] .= 0.23 .* Ct[:,:,ID_Turbines].^-0.25 .* TI_0_vec[:,:,ID_Turbines].^0.17
    a[:,:,ID_Turbines]       .= 0.93 .* Ct[:,:,ID_Turbines].^-0.75 .* TI_0_vec[:,:,ID_Turbines].^0.17
    b[:,:,ID_Turbines]       .= 0.42 .* Ct[:,:,ID_Turbines].^0.6   .* TI_0_vec[:,:,ID_Turbines].^0.2
    c[:,:,ID_Turbines]       .= 0.15 .* Ct[:,:,ID_Turbines].^-0.25 .* TI_0_vec[:,:,ID_Turbines].^-0.7
    d[:,:,ID_Turbines]       .= 2.3  .* Ct[:,:,ID_Turbines].^1.2   .* TI_0_vec[:,:,ID_Turbines].^0.1
    e[:,:,ID_Turbines]       .= 1.0                                   .* TI_0_vec[:,:,ID_Turbines].^0.1
    f[:,:,ID_Turbines]       .= 0.7  .* Ct[:,:,ID_Turbines].^-3.2  .* TI_0_vec[:,:,ID_Turbines].^-0.45

    println(sum(ID_Turbines))

    return k, epsilon, a, b, c, d, e, f
end#ComputeEmpiricalVars
    
function comp_ConvectionVelocity(tmp, u_c_vec, Ct, D, sigma, u_0_vec)
# Compute local convection velocity of each turbine
    # Term within the sqrt. expression
    tmp = (1 .- ((Ct.*D^2)./(8 .* sigma[:,1:1,:].^2)))
    tmp[tmp .< 0] .= 0          #Correction (-Inf appears for those points which should not be computed, since they are upstream. Additionally negative numbers appear due to incompatibility between Ishihara-Qian & convection velocity derivation for the near wake.)
    tmp[isnan.(tmp)] .= 0
    u_c_vec .= (0.5 .+ 0.5 .* sqrt.(tmp)) .* u_0_vec
    #u_c_vec[tmp .== 0] .= 0     #Second correction, to prevent bugs when used for superposition
    return u_c_vec
end#comp_ConvectionVelocity
    
function Meandering_Correction(sigma_m::Array{Float64}, psi::Array{Float64}, Lambda::Array{Float64}, TI_0_vec::Array{Float64}, u_0_vec::Array{Float64}, ZCoordinates::Array{Float64}, XCoordinates::Array{Float64}, u_c_vec::Array{Float64}, Call_for_Uc::Bool)
# Compute the meandering correction according to the Braunbehrens & Segalini model (2019) 
    if Call_for_Uc == false #Call only during rotor point computation
    #Compute fluctuation intensity
    psi     .= 0.7.*TI_0_vec.*u_0_vec
    end

    #Compute integral length scale of representative eddy   
    Lambda  .= (0.4 .* ZCoordinates) ./ psi
    #Compute corrected wake width 
    sigma_m .= sqrt.((2 .* psi .* Lambda.^2) .* ((XCoordinates./(u_c_vec.*Lambda)) .+ exp.(-XCoordinates./(u_c_vec.*Lambda)) .- 1))  
    return psi, Lambda, sigma_m          
end#Meandering_Correction
    
#####################################################
end