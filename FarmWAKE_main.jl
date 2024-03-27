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


x=1

