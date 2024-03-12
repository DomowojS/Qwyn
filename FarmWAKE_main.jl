#=  This Script excexutes the FarmWAKE tool. 
    Inputs are taken from 01_Inputs.
    All function blocks can be found in 02_Modules
=#

#Excecute PKG script, add/ uodate all neccecary packages can be turned of after first run
include("02_Modules/PKG_Manager.jl")
#Run PKG_Manager, when in doubt if you have all necesarry packages installed.
#PKG_Manager()

#Generate all Input file Structs & Initialise Matrices
include("02_Modules/Preprocessing.jl")

generateWF("WF", "01_Input")

