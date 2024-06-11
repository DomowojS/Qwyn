# Creates a dictionary of Input data
using JLD2 #Needed for loading turbine & meteorological data
using LinearAlgebra
userdata = OrderedDict{String, Any}(
    ##########      (1) Wind farm data         ######################
    # Name of the wind Farm
    "name" => "HR1",
    # Wind Farm Data
    "N" => 80, #Number of turbines  
    "x_vec" => [0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 0, 7, 14, 21, 28, 35, 42, 49, 56, 63], #LinRange(0,21,80),   #X-Coordinates [0, 7, 14, 21, 28, 35],
    "y_vec" => [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49], # zeros(80,),         #Y-Coordinates [0, 0, 0, 0, 0, 0],

    ##########      (2) Turbine data           ######################
    "Yaw" => 270 .+ zeros(80,), # Yaw angle of the turbines (In geographical DEG)
    #Turbine Type
    "Turbine_Type" => "VestasV80",  #Turbine Type. One type for the whole wind farm. 
                                    #Possible Inputs: "VestasV80", "NREL_5MW", "DTU_10MW", "IEA_15MW"

    ##########      (3) Atmospheric data       ######################
    #Use either 3.1 for single computation OR 3.2 for AEP computation
    # (3.1) Single computatiosn 
    #       This section is only used for single case computation    
    "u_ambient" => 14,     # [m/s] Ambient wind speed
    "alpha"     => 270,     # [°] Geographical direction of the wind. -> N == 0°
    "TI_a"      => 0.07,     # [-] Ambient turbulence intensity in [-]
    "z_Surf"    => 0.005,     # [-] Surface roughness of the modelled case *for offshore conditions z_Surf should equal between 0.0001 (calm see) and 0.01 (high waves)
    "z_r"       => 70.0,    # [m] Height the average wind speed "u_ambient" was measured. If not known, choose z_r = 10

    # (3.2) AEP computation 
    #       This section is only used for AEP computation      
    "Wind_rose" => 3.14159999, # Get wind rose as specified in "04_Ambient_data"

    ##########      (4) Computational setting  ######################
    "SimpleComp"    => true,    # For the computation of one case 
    "AEPComp"       => false,   # For the estimation of the farms AEP
    ## Advanced settings:
    #Superposition Method
    "Superpos"  => "Linear_Rotorbased", #Superposition method for velocity deficits. Choose between linear rotorbased summation & momentum conserving approach. 
                                           #Possible inputs: "Linear_Rotorbased", "Momentum_Conserving"

    ##########      (5) Numerical parameters   ######################
    "Dimensions"            => "3D",        #Choose dimensions resolution. 1) Three dimensional space or 2) two dimansional plane at Hub height.
                                            #Possible inputs: "3D", "2D"
    "Rotor_Discretization"  => "gridded",   #Specifies the rotor descritization technique. Current choices: 1) Evenly distributed grid (slow with small error), 2) Fibonacci-Latice distributed points (quicker). 
                                            #Possible inputs: "gridded", "fibonacci"
    "Rotor_Res"             => 100,         #Number of points used to represent the rotor. Reccomendation: 100 for "griddeed" & XX for "fibonacci".

    ##########      (6) Graphical output       ######################
    "z" => 100,    
    "Z_Max"     => 110,  #Maximum height
    "Z_Min"     => 30,   #Minimum height

    ##########   Literature Input              ######################
    #= These Numbers are placeholders. They overwritten by data from literature.
    You can provide this as .jld2 file for turbine & atmospheric data.
        1) Turbine data     --> provide .jld2 file providing turbine data in folder "04_Turbine_Data"
        2) Atmospheric data --> provide .jld2 file providing atmospheric data in folder "05_Atmospheric_Data"
    A documentation on how to provide the data can be found in the corresponding folders.
    =#
    #Turbine placeholders:
    "D"         => 0.0,        # Turbine diameter 
    "H"         => 0.0,        # Hub height in [m]
    "P_Input"   => zeros(1,1), # Power coefficient - defined as .txt in "03_Turbine_Data"
    "Ct_Input"  => zeros(1,1), # Thrust coefficient - defined as .txt in "03_Turbine_data"
    #Atmospheric data placeholders:
    "u_ambient_zprofile" => zeros(1,1,1), # [m/s] height profile of the wind as vector of z coordinates resulting from amount of rotor resolution points 
)
