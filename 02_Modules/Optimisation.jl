#=
Optimisation module
=#

module Optimisation
include("Input_Processing.jl")
include("Initialisation_Module.jl")
include("SimpleComputation.jl")
using .SimpleComputation, .Input_Processing, .Initialisation_Module, NLopt

export SetAndRunOptimiser, AEP4Layout_Optim

function SetAndRunOptimiser(WindFarm, speeds, angles, frequencies, Powers_AllCases)
#Function to prepare and set and run the optimisation model
# 1) Variables are set & starting values are defined
# 2) Constraint functions are applied
# 3) Handle to target function is passed
# 4) Optimisation problem is run
# 5) Optimised parameters are returned

    n = WindFarm.N      # Counter for amount of turbines
    min_distance = 3.7  # Minimum distance between turbines --> To be exported to config file

    # Optimiser settings

    opt = NLopt.Opt(:LN_COBYLA, WindFarm.N*2)  # 2*Turbine number design variables (x & y coordinates)        
    #NLopt.lower_bounds!(opt, vcat(fill(minimum(WindFarm.x_vec), n), fill(minimum(WindFarm.y_vec), n)))      # Lower boundaries (Box condition)
    #NLopt.upper_bounds!(opt, vcat(fill(maximum(WindFarm.x_vec), n), fill(maximum(WindFarm.y_vec), n)))      # Upper boundaries (Box condition)

    # Only For Horns Rev!
    NLopt.lower_bounds!(opt, vcat(fill(-Inf, n), fill(minimum(WindFarm.y_vec), n)))
    NLopt.upper_bounds!(opt, vcat(fill(Inf, n), fill(maximum(WindFarm.y_vec), n)))
    
    
    #Tolerance setting
    NLopt.xtol_rel!(opt, 1e-4)  
    
    #Set Objective function
    NLopt.min_objective!(opt, (x, grad) -> objective_func(x, grad, WindFarm, speeds, angles, frequencies, Powers_AllCases))
    
    # Add pairwise distance constraints
    for i in 1:WindFarm.N
        for j in (i+1):WindFarm.N
            NLopt.inequality_constraint!(opt, 
                (x, g) -> distance_constraint(x, g, i, j, min_distance), 
                1e-8)
        end
    end

    # Only for Horns Rev1
    # Add slope constraints for x coordinates
    for i in 1:n
        # Lower bound slope constraint: x >= -6/48.63 * y
        NLopt.inequality_constraint!(opt,
            (x, g) -> slope_constraint_lower(x, g, i, n),
            1e-8)
        
        # Upper bound slope constraint: x <= (-6/48.63 * y) + 69
        NLopt.inequality_constraint!(opt,
            (x, g) -> slope_constraint_upper(x, g, i, n),
            1e-8)
    end

    min_f, min_x, ret = NLopt.optimize(opt, vcat(WindFarm.x_vec, WindFarm.y_vec))
    num_evals = NLopt.numevals(opt)


    # Extract results
    OptimalAEP = -min_f
    x_vec_optimal = min_x[1:n]
    y_vec_optimal = min_x[n+1:2n]
    

    return OptimalAEP, x_vec_optimal, y_vec_optimal, ret, num_evals

end#SetAndRunOptimiser

function AEP4Layout_Optim(WindFarm, speeds, angles, frequencies, Powers_AllCases::Vector{Float64})
#=  Target function for AEP Optimisation.
    Computes AEP without comments by including the optimisers suggested parameters. 

    ToDo:   Include Flag option to replace the "right" parameters (YawAngle/Layout etc.)
            so this function becomes reusable for all cases
=#
    
    ######## START LOOP OVER all wind speeds HERE ############################
        
        for case in 1:length(frequencies)

            # Assign right wind speed & angle for this iteration
            WindFarm.u_ambient = speeds[case] #speed
            WindFarm.alpha     = angles[case] #angle

            #Initialise all arrays & matrices needed for the computation.
            WindFarm, CS = initCompArrays(WindFarm) #Initialises mutable struct "CA" which contains necessary 
                                                    #computation arrays & computes coordinates acc. to user Input.

            LoadTurbineDATA!(WindFarm, CS)          #Update Input & computation structs with provided power & 
                                                    #thrust curves
            
            WindFarm.u_ambient_zprofile, CS = LoadAtmosphericData(WindFarm,CS)      #Update Input & computation structs with atmospheric data &
                                                                                    #(wind shear profile, wind rose etc.)
            
            FindStreamwiseOrder!(WindFarm, CS)



            # Iterating over turbine rows
            for CS.i=1:maximum(CS.CompOrder)
                FindComputationRegion!(WindFarm, CS)
                
                if any(CS.ID_Turbines != 0) # Check if any turbine's wake has to be computed

                Single_Wake_Computation!(WindFarm, CS)  #Compute single wake effect
                
                Superposition!(WindFarm, CS, WindFarm.u_ambient_zprofile)            #Compute mixed wake
                
                getTurbineInflow!(WindFarm, CS)         #Evaluate new inflow data
                
                getNewThrustandPower!(WindFarm, CS)     #Evaluate new operation properties

                end
            end
                
            getTotalPower!(CS)  #Compute total power of the wind farm 

            Powers_AllCases[case]=CS.TotalPower;   #Save computed power for this case

        end
            
            TotalAEP=sum(Powers_AllCases.*frequencies) * 8760;               #Compute total AEP by multyplying power with frequencies
     
        return TotalAEP

end#AEP4Layout_Optim


function objective_func(x::Vector, grad::Vector, WindFarm, speeds, angles, frequencies, Powers_AllCases)
#=  Defines Objective function for Optimization package
    Pure definition/ assignment block, no further functionality
=#
    #Overwrite x_vec and y_vec vectors with optimiser's variables
    n = length(x) ÷ 2
    WindFarm.x_vec = x[1:n]
    WindFarm.y_vec = x[n+1:2n]


    return -AEP4Layout_Optim(WindFarm, speeds, angles, frequencies, Powers_AllCases)
end#objective_func!

function distance_constraint(x::Vector, grad::Vector, i::Int, j::Int, max_dist::Float64)
#=  Defines minimum distance constraints between turbines.
    Acts like "Safety circle around turbines, so that they can't be placed
    in each others nearwake.
=#

    dist = get_machine_distance(x, i, j)
    return max_dist - dist

end#distance_constraints!


function get_machine_distance(x::Vector, i::Int, j::Int)
# Helper function to compute distance between two machines
    N = length(x) ÷ 2
    # Get coordinates for machine i
    xi, yi = x[i], x[i + N]
    # Get coordinates for machine j
    xj, yj = x[j], x[j + N]
    # Compute Euclidean distance
    return sqrt((xi - xj)^2 + (yi - yj)^2)
end#get_machine_distance


function slope_constraint_lower(x, g, i, n)
    # x[i] >= -6/48.63 * x[i+n]
    # Rearrange to: (-6/48.63 * y) - x <= 0
    return (-6/48.63 * x[i+n]) - x[i]  # Returns <= 0 when constraint is satisfied
end

function slope_constraint_upper(x, g, i, n)
    # x[i] <= (-6/48.63 * x[i+n]) + 69
    # Already in form: x - (-6/48.63 * y) - 69 <= 0
    return x[i] - (-6/48.63 * x[i+n]) - 69  # Returns <= 0 when constraint is satisfied
end

end