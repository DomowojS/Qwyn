#=  This Script excexutes the FarmWAKE tool. 
    Inputs are taken from 01_Inputs.
    All function blocks can be found in 02_Modules
=#

#Excecute PKG script, adding all neccecary packages
include("02_Modules/PKG_Manager.jl")
PKG_Manager()
#Generate all Input file Structs & Initialise Matrices
include("02_Modules/generateWF.jl")


