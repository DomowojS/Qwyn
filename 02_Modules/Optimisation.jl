#=
Optimisation module
=#

module Optimisation
include("Input_Processing.jl")
include("Initialisation_Module.jl")
include("SimpleComputation.jl")
using .SimpleComputation, .Input_Processing, .Initialisation_Module, Optimisation

export SetAndRunOptimiser

function SetAndRunOptimiser(WindFarm, frequencies)
#Function to prepare and set and run the optimisation model
# 1) Variables are set & starting values are defined
# 2) Constraint functions are applied
# 3) Handle to target function is passed
# 4) Optimisation problem is run
# 5) Optimised parameters are returned

    n = WindFarm.N      # Counter for amount of turbines
    min_distance = 3.7  # Minimum distance between turbines --> To be exported to config file
    
    # Create model
    Optimise_Layout = Model(Ipopt.Optimizer)

    # Define variables
    @variable(Optimise_Layout, x[1:n])
    @variable(Optimise_Layout, y[1:n])

    # Set start values
    for i in 1:n
        set_start_value(x[i], WindFarm.x_vec[i])
        set_start_value(y[i], WindFarm.y_vec[i])
    end

    # Set any bounds if needed --> To be exported to config file
    for i in 1:n
        set_lower_bound(y[i], minimum(WindFarm.y_vec))
        set_upper_bound(y[i], maximum(WindFarm.y_vec))
    end

    # Add minimum distance constraints between all turbines
    for i in 1:n
        for j in (i+1):n
            @NLconstraint(Optimise_Layout, 
                sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2) >= min_distance
            )
        end
    end

    # Add x constraint
    for i in 1:n
        set_lower_bound(x[i], minimum(WindFarm.x_vec))
        set_upper_bound(x[i], maximum(WindFarm.x_vec))
        #@constraint(Optimise_Layout, x[i] == x_constraint(y[i]))     # Try later with HR1 slopes!
    end

    # Define ojective function & optimisation function
    @objective(Optimise_Layout, Max, AEP4Layout_Optim(x, y, WindFarm, frequencies))

    # Solve
    optimize!(Optimise_Layout)

    # Get results
    x_vec_optimal = value.(x)
    y_vec_optimal = value.(y)
    OptimalAEP = objective_value(Optimise_Layout)

    return OptimalAEP, x_vec_optimal, y_vec_optimal

end#SetAndRunOptimiser

function AEP4Layout_Optim(x, y, WindFarm, frequencies)
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



end