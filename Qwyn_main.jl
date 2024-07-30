include("02_Modules/PKG_Manager.jl")            #pck manager (for own overiew -> !!to be deleted!! 
include("02_Modules/Input_Processing.jl")       #Module to process input data
include("02_Modules/Initialisation_Module.jl")  #Module for array initialisation/ space preallocation
include("02_Modules/SimpleComputation.jl")      #Module for simple computation.
include("02_Modules/Postprocessing.jl")
using .Input_Processing, .Initialisation_Module, .SimpleComputation, .Postprocessing, TickTock, MAT,Plots

function Qwyn_Simple(u_ambient::Real, alpha::Real, TI_a::Real)
#=  This Script excexutes Qwyn. 
    Inputs are taken from 01_Inputs.
    All function blocks can be found in 02_Modules
=#
    tick()
    WF=generateWF("WF", "01_Inputs", u_ambient, alpha, TI_a) #Generate array consisting of all input data for each Input 
                                                             #file in "01_Inpts"
    
    # Iterate over defined cases. Compute one by one (according to settings)
    for WindFarm in WF
        println("###########################")
        println("Computing: ", WindFarm.name)   #Terminal output for which input file is being processed
        println("###########################")

        #Initialise all arrays & matrices needed for the computation.
        println("Initialising arrays..")
        WindFarm, CS = initCompArrays(WindFarm) #Initialises mutable struct "CA" which contains necessary 
                                                #computation arrays & computes coordinates acc. to user Input.

        println("..loading turbine data..")
        LoadTurbineDATA!(WindFarm, CS)          #Update Input & computation structs with provided power & 
                                                #thrust curves
        
        println("..loading atmospheric data..")
        LoadAtmosphericData!(WindFarm,CS)       #Update Input & computation structs with atmospheric data &
                                                #(wind shear profile, wind rose etc.)
        
        println("..finding computation order..")
        FindStreamwiseOrder!(WindFarm, CS)

        tock()
        tick()

        # Iterating over turbine rows
        println("Starting iterative computation...")
        for CS.i=1:maximum(CS.CompOrder)
            println("Iteration ", CS.i)

            FindComputationRegion!(WindFarm, CS)
            
            if any(CS.ID_Turbines != 0) # Check if any turbine's wake has to be computed

            Single_Wake_Computation!(WindFarm, CS)  #Compute single wake effect
            
            Superposition!(WindFarm, CS)            #Compute mixed wake
            
            getTurbineInflow!(WindFarm, CS)         #Evaluate new inflow data
            
            getNewThrustandPower!(WindFarm, CS)     #Evaluate new operation properties

            end
        end
            
        getTotalPower!(CS)  #Compute total power of the wind farm 
        println("...computation finished!")
        tock()

        println("Plots on the way...")
        SimplePlots(WindFarm, CS)
        println("...finished!")
    
    return CS
    #TMP=reshape(CS.P_vec[[4, 12, 20, 28, 36, 44, 52, 60]]./CS.P_vec[4], 8) #270
    #TMP=reshape(CS.P_vec[[5, 12, 19, 26, 33]]./CS.P_vec[4], 5)              #222
    TMP=reshape(CS.P_vec[[4, 13, 22, 31, 40]]./CS.P_vec[4], 5)             #312
    Base.print_matrix(stdout, TMP)
    end
    return WF;
end#Qwyn_Simple

function Qwyn_AEP()
#= -> Logic: To be done -> Using different functions for different tasks.

    1) Qwyn_Simple
    2) Qwyn_AEP

    ( 3) Qwyn_Optimise_AEP
    4) Qwyn_Optimise_Yaw ) => Could be written by a separate user 

    Idea: Possible to define a functon which can go without input (default) or with input -> when called by another function.
    Meaning: If the user wants to compute an AEP. The AEP function simply calls the "Qwyn_Simple" function several times.

    Followup Idea: 
    Define a second "Qwyn_Simple". For example "Qwyn_Simple_internal_use" --> can be called by other function such as "Qwyn_AEP". 
    Here, the structs WindFarm & CS can be passed directly. without having to process the input data again and initalise matrices.

    Meaning for Input file:
    The Input file will, hence, be reduced to physical & numerical settings.
=#
end

function Qwyn_Optimiser()
# Optimisation functions will be added here.

end#Qwyn_AEP()

#...more function to come?























                                #= Temporary stuff -> for saving to MAT and plotting in Matlab!    
                                    # Convert the struct to a dictionary
                                    CS.Interp_Ct=0
                                    CS.Interp_P=0
                                struct_dict = Dict{String, Any}(string.(propertynames(CS)) .=> getfield.(Ref(CS), propertynames(CS)))
                                struct_dict2 = Dict{String, Any}(string.(propertynames(WindFarm)) .=> getfield.(Ref(WindFarm), propertynames(WindFarm)))

                                # Specify the filename for the .mat file
                                filename = "99_PlotWMATLAB/WindFarmCS_NewCoordinateSys.mat"
                                filename2 = "99_PlotWMATLAB/WindFarmWF_NewCoordinateSys.mat"
                                # Save the struct to the .mat file
                                matwrite(filename, struct_dict)
                                matwrite(filename2, struct_dict2)
                                =#