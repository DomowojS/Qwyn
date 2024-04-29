# Creates a dictionary of Input data
using JLD2 #Needed for loading turbine & meteorological data
userdata = OrderedDict{String, Any}(
    ##########      (1) Wind farm data         ######################
    # Name of the wind Farm
    "name" => "HR1",
    # Wind Farm Data
    "N" => 2, #Number of turbines  
    "x_vec" => [0.0, 7.0],   #X-Coordinates
    "y_vec" => [0.0, 0.0],     #Y-Coordinates

    ##########      (2) Turbine data           ######################
    "Yaw" => [270, 270], # Yaw angle of the turbines (In geographical DEG)
    #Turbine Type
    "VestasV80" => true,
    "NREL_5MW"  => false,

    # Data from literature -> provide as .jld2 file in folder "04_Turbine_Data"
    # All following turbine data are placeholders & will be loaded according to choice of turbine
    "D" => 0.0, # Turbine diameter 
    "H" => 0.0, # Hub height in [m]
    "P_Input"  => zeros(2,2), # Power coefficient - defined as .txt in "03_Turbine_Data"
    "Ct_Input" => zeros(2,2), # Thrust coefficient - defined as .txt in "03_Turbine_data"

    ##########      (3) Atmospheric data       ######################
    #Use either 3.1 for single computation OR 3.2 for AEP computation
    # (3.1) Single computatiosn 
    #       This section is only used for single case computation    
    "u_ambient" => 5.0, # Ambient wind speed in [m/s]
    "alpha" => 270,    # Geographical direction of the wind speed in [°]. -> N == 0°
    "TI_a" => 0.07,    # Ambient turbulence intensity in [-]
    # (3.2) AEP computation 
    #       This section is only used for AEP computation      
    "Wind_rose" => 3.14159999, # Get wind rose as specified in "04_Ambient_data"

    ##########      (4) Computational setting  ######################
    "SimpleComp" => true,
    "AEPComp" => false,
    "Optim" => false,

    ##########      (5) Numerical parameters   ######################
    "RotorRes" => 1000, #Number of points used to distrectisize the turbine's rotors

    ##########      (6) Graphical output       ######################
    "z" => 100

)
