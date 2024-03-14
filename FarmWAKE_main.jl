#=  This Script excexutes the FarmWAKE tool. 
    Inputs are taken from 01_Inputs.
    All function blocks can be found in 02_Modules
=#
#Excecute PKG script, add/ uodate all neccecary packages can be turned of after first run
include("02_Modules/PKG_Manager.jl")

#Generate all Input file Structs & Initialise Matrices
include("02_Modules/Preprocessing.jl")  #Include module with relevant functions
import .Preprocessing: generateWF                #Make them available in this script.

WF=generateWF("WF", "01_Inputs")        #



