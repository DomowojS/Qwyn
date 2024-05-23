# Creates a dictionary of Input data
using JLD2 #Needed for loading turbine & meteorological data
using LinearAlgebra
userdata = OrderedDict{String, Any}(
    ##########      (1) Wind farm data         ######################
    # Name of the wind Farm
    "name" => "HR1",
    # Wind Farm Data
    "N" => 80, #Number of turbines  
    "x_vec" => LinRange(0,21,80), #[0, 7, 14],   #X-Coordinates
    "y_vec" => zeros(80,),     #Y-Coordinates

    ##########      (2) Turbine data           ######################
    "Yaw" => 270 .+ zeros(80,), # Yaw angle of the turbines (In geographical DEG)
    #Turbine Type
    "VestasV80" => true,
    "NREL_5MW"  => false,

    ##########      (3) Atmospheric data       ######################
    #Use either 3.1 for single computation OR 3.2 for AEP computation
    # (3.1) Single computatiosn 
    #       This section is only used for single case computation    
    "u_ambient" => 10.0,     # [m/s] Ambient wind speed
    "alpha"     => 270,     # [°] Geographical direction of the wind speed. -> N == 0°
    "TI_a"      => 0.1,     # [-] Ambient turbulence intensity in [-]
    "z_Surf"    => 0.005,   # [-] Surface roughness of the modelled case *for offshore conditions z_Surf should equal between 0.0001 (calm see) and 0.01 (high waves)
    "z_r"       => 70.0,    # [m] Height the average wind speed "u_ambient" was measured. If not known, choose z_u = 10

    # (3.2) AEP computation 
    #       This section is only used for AEP computation      
    "Wind_rose" => 3.14159999, # Get wind rose as specified in "04_Ambient_data"

    ##########      (4) Computational setting  ######################
    "SimpleComp"        => true,    # For the computation of one case 
    "AEPComp"           => false,   # For the estimation of the farms AEP
    ## Advanced settings:
    #Superposition Method
    "Linear_Rotorbased"     => true,   # Superposition using linear rotorbased summation for velocity deficit
    "Momentum_Conserving"   => false,  # Superposition using momentum conserving approach for velocity deficit

    ##########      (5) Numerical parameters   ######################
    "Y_Res"     => 100,     #Number of spanwise points used to distrectisize the turbine's rotors
    "Z_Res"     => 10,      #Number of height points to descritisise the rooms (number of height levels computed)

    ##########      (6) Graphical output       ######################
    "z" => 100,
    "Z_Max"     => 1.5*80,  #Maximum height
    "Z_Min"     => 0,       #Minimum height
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
    "u_ambient_zprofile" => zeros(1,1,1,1), # [m/s] height profile of the wind as vector of z coordinates resulting from amount of rotor resolution points 
)
