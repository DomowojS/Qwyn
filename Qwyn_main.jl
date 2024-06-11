using TickTock, MAT
#Excecute PKG script, add/ update all neccecary packages can be turned of after first run
include("02_Modules/PKG_Manager.jl")
#Read all input files & create struct array of the data
include("02_Modules/Input_Processing.jl")  #Precompile module with relevant functions
#Compute / Optimise as requested in each input file 
include("02_Modules/Initialisation_Module.jl") #Include module which contains initialisation of computational arrays -> Needed for every operation.
include("02_Modules/SimpleComputation.jl")
using .SimpleComputation
using .Input_Processing
using .Initialisation_Module

function Qwyn()
#=  This Script excexutes Qwyn. 
    Inputs are taken from 01_Inputs.
    All function blocks can be found in 02_Modules
=#
    tick()
    WF=generateWF("WF", "01_Inputs")           #Generate array consisting of all input data for each Input file in "01_Inptus"
    
    # Iterate over defined cases. Compute one by one (according to settings)
    for WindFarm in WF
        println("###########################")
        println("Computing: ", WindFarm.name)   #Terminal output for which input file is being processed
        println("###########################")

        #Initialise all arrays & matrices needed for the computation.
        println("Initialising arrays..")
        WindFarm, CS = initCompArrays(WindFarm) #Initialises mutable struct "CA" which contains necessary computation arrays & computes coordinates acc. to user Input.
        
        println("..loading turbine data..")
        LoadTurbineDATA!(WindFarm, CS)          #Update Input & computation structs with provided power & thrust curves
        println("..loading atmospheric data..")
        LoadAtmosphericData!(WindFarm,CS)       #Update Input & computation structs with atmospheric data (wind shear profile, wind rose etc.)
        tock()
        tick()
        # Iterating over turbine rows
        i=0
        println("Starting iterative computation...")
        while CS.zeta > 10^-3
            i=i+1
            println("Iteration ", i)
            #Compute single wake effect
            Ishihara_WakeModel!(WindFarm, CS)  #Compute wakes of single turbines
            #Compute mixed wake
            Superposition!(WindFarm, CS)
            #Evaluate new inflow data
            getTurbineInflow!(WindFarm, CS) 
            #Evaluate new operation properties
            getNewThrustandPower!(WindFarm, CS)
        end
            #Compute total power of the wind farm 
            getTotalPower!(CS)
        println("...computation finished!")
        tock()



    ## Temporary stuff     
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

    global CS
    end

    return WF, CS;
end#Qwyn


