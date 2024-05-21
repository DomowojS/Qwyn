#=  This Script excexutes the Qwyn. 
    Inputs are taken from 01_Inputs.
    All function blocks can be found in 02_Modules
=#
#Excecute PKG script, add/ update all neccecary packages can be turned of after first run
include("02_Modules/PKG_Manager.jl")
using TickTock
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
    println("Computing: ", WindFarm.name) #Terminal output for which input file is being processed
    println("###########################")
    tick()
    #Initialise all arrays & matrices needed for the computation.
    WindFarm, CS = initCompArrays(WindFarm)         #Initialises mutable struct "CA" which contains necessary computation arrays & computes coordinates acc. to user Input.
    LoadTurbineDATA!(WindFarm, CS)    #Update Input & computation structs with provided power & thrust curves
    LoadAtmosphericData!(WindFarm,CS) #Update Input & computation structs with atmospheric data (wind shear profile, wind rose etc.)
    tock()
    tick()
    #Compute single wake effect
    Ishihara_WakeModel!(WindFarm, CS)  #Compute wakes of single turbines
    tock()
    global CD = CS
end
println("This was a Test")

