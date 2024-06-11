#=  This Script excexutes the Qwyn. 
    Inputs are taken from 01_Inputs.
    All function blocks can be found in 02_Modules
=#
#Excecute PKG script, add/ update all neccecary packages can be turned of after first run
include("02_Modules/PKG_Manager.jl")
using TickTock, MAT
tick()
#Read all input files & create struct array of the data
include("02_Modules/Input_Processing.jl")  #Precompile module with relevant functions
using .Input_Processing                    #Make all contents of the module available in this script (Import only imports module name)

WF=generateWF("WF", "01_Inputs")           #Generate array consisting of all input data for each Input file in "01_Inptus"

#Compute / Optimise as requested in each input file 
include("02_Modules/Initialisation_Module.jl") #Include module which contains initialisation of computational arrays -> Needed for every operation.
using .Initialisation_Module

if any(setting.SimpleComp for setting in WF) #Include only modules which are needed (according to serttings)
    include("02_Modules/SimpleComputation.jl")
    using .SimpleComputation
elseif any(setting.AEPComp for setting in WF)
    include("02_Modules/SimpleComputation.jl")
    include("02_Modules/AEPComputation.jl")
    using .SimpleComputation
    using .AEPComputation
elseif any(setting.Optimisation for setting in WF)
    include("02_Modules/SimpleComputation.jl")
    include("02_Modules/AEPComputation.jl")
    include("02_Modules/Optimisation.jl")
    using .SimpleComputation
    using .AEPComputation
    using .Optimisation
else
    println("No computation request in any Input file. Task finalised.")
end
    
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
    println("...computation finished!")
    tock()
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

global CD = CS

end

