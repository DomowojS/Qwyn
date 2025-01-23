#=
Optimisation module
=#

module Optimisation
include("Input_Processing.jl")
include("Initialisation_Module.jl")
include("SimpleComputation.jl")
using .SimpleComputation, .Input_Processing, .Initialisation_Module, Optimization, Ipopt, FiniteDiff, GalacticOptim

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

    # Initial guess taken from User Inputs
    x0 = vcat(WindFarm.x_vec, WindFarm.y_vec)

    # Bounds
    lb = vcat(
        #fill(minimum(WindFarm.x_vec)*4, n), 
        #fill(minimum(WindFarm.y_vec)*4, n)
        fill(0, n), 
        fill(0, n)
    )
    ub = vcat(
        fill(20, n), 
        fill(5, n)
    )

    # Number of distance constraints
    num_constraints = div(n * (n-1), 2)

    # Optimization problem setup
    optf = OptimizationFunction(
        (F,x) -> objective_func!(F, x, WindFarm, speeds, angles, frequencies, Powers_AllCases),
        GalacticOptim.AutoFiniteDiff(), 
        cons = (c, x) -> distance_constraints!(c, x, min_distance)
    )

    optprob = OptimizationProblem(
        optf, 
        x0, 
        lb = lb, 
        ub = ub
        #lcons = zeros(num_constraints),     # Lower bounds for constraints (c[k] >= 0)
        #ucons = fill(Inf, num_constraints)  # No upper bounds (c[k] unbounded above)
    )

    # Solve
    sol = solve(optprob, GalacticOptim.AugmentedLagrangian())

    # Extract results
    x_vec_optimal = sol.u[1:n]
    y_vec_optimal = sol.u[n+1:2n]
    OptimalAEP = -sol.minimum

    return OptimalAEP, x_vec_optimal, y_vec_optimal

end#SetAndRunOptimiser

function AEP4Layout_Optim(x, y, WindFarm, speeds, angles, frequencies, Powers_AllCases::Vector{Float64})
#=  Target function for AEP Optimisation.
    Computes AEP without comments by including the optimisers suggested parameters. 

    ToDo:   Include Flag option to replace the "right" parameters (YawAngle/Layout etc.)
            so this function becomes reusable for all cases
=#

    #Overwrite x_vec and y_vec vectors with optimiser's variables
    tempWindFarm = deepcopy(WindFarm)
    tempWindFarm.x_vec = x
    tempWindFarm.y_vec = y
    
    ######## START LOOP OVER all wind speeds HERE ############################
        
        for case in 1:length(frequencies)

            # Assign right wind speed & angle for this iteration
            tempWindFarm.u_ambient = speeds[case] #speed
            tempWindFarm.alpha     = angles[case] #angle

            #Initialise all arrays & matrices needed for the computation.
            tempWindFarm, CS = initCompArrays(tempWindFarm) #Initialises mutable struct "CA" which contains necessary 
                                                    #computation arrays & computes coordinates acc. to user Input.

            LoadTurbineDATA!(tempWindFarm, CS)          #Update Input & computation structs with provided power & 
                                                    #thrust curves
            
            tempWindFarm.u_ambient_zprofile, CS = LoadAtmosphericData(tempWindFarm,CS)      #Update Input & computation structs with atmospheric data &
                                                                                    #(wind shear profile, wind rose etc.)
            
            FindStreamwiseOrder!(tempWindFarm, CS)



            # Iterating over turbine rows
            for CS.i=1:maximum(CS.CompOrder)
                FindComputationRegion!(tempWindFarm, CS)
                
                if any(CS.ID_Turbines != 0) # Check if any turbine's wake has to be computed

                Single_Wake_Computation!(tempWindFarm, CS)  #Compute single wake effect
                
                Superposition!(tempWindFarm, CS, tempWindFarm.u_ambient_zprofile)            #Compute mixed wake
                
                getTurbineInflow!(tempWindFarm, CS)         #Evaluate new inflow data
                
                getNewThrustandPower!(tempWindFarm, CS)     #Evaluate new operation properties

                end
            end
                
            getTotalPower!(CS)  #Compute total power of the wind farm 

            Powers_AllCases[case]=CS.TotalPower;   #Save computed power for this case

        end
            
            TotalAEP=sum(Powers_AllCases.*frequencies) * 8760;               #Compute total AEP by multyplying power with frequencies
     
        return TotalAEP

end#AEP4Layout_Optim


function objective_func!(F, x, WindFarm, speeds, angles, frequencies, Powers_AllCases)
#=  Defines Objective function for Optimization package
    Pure definition/ assignment block, no further functionality
=#
    n = length(x) ÷ 2
    x_vec = x[1:n]
    y_vec = x[n+1:2n]
    F[1] = -AEP4Layout_Optim(x_vec, y_vec, WindFarm, speeds, angles, frequencies, Powers_AllCases)
end#objective_func!

function distance_constraints!(c, x, min_distance)
#=  Defines minimum distance constraints between turbines.
    Acts like "Safety circle around turbines, so that they can't be placed
    in each others nearwake.
=#
    n = length(x) ÷ 2
    x_vec = x[1:n]
    y_vec = x[n+1:2n]
    k = 1
    for i in 1:n
        for j in (i+1):n
            c[k] = sqrt((x_vec[i] - x_vec[j])^2 + (y_vec[i] - y_vec[j])^2) - min_distance
            k += 1
        end
    end
    c[n+1:end]=c[1:n]
end#distance_constraints!

end