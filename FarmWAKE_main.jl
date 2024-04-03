#=  This Script excexutes the FarmWAKE tool. 
    Inputs are taken from 01_Inputs.
    All function blocks can be found in 02_Modules
=#
#Excecute PKG script, add/ update all neccecary packages can be turned of after first run
include("02_Modules/PKG_Manager.jl")

#Generate all Input file Structs & Initialise Matrices
include("02_Modules/Preprocessing.jl")  #Precompile module with relevant functions
using .Preprocessing                    #Make all contents of the module available in this script (Import only imports module name)

WF=generateWF("WF", "01_Inputs")        #Generate array consisting of all input data for each Input file in "01_Inptus"

#Compute / Optimise as requested in each input file 
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


println(WindFarm.name)

end


x=1

