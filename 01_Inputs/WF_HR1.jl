# Creates a dictionary of Input data
userdata = OrderedDict{String, Any}(
    ##########      (1) Wind farm data         ######################
    # Name of the wind Farm
    "name" => "HR1",
    # Wind Farm Data
    "N" => 3, #Number of turbines  
    "x_vec" => [0.0, 7.0, 14.0],   #X-Coordinates
    "y_vec" => [0.0, 0.0, 7.0],     #Y-Coordinates

    ##########      (2) Turbine data           ######################
    "D" => 80.0,             # Turbine diameter in [m]
    "H" => 70.0,             # Hub height in [m]
    "Cp" => [1.0, 2.0],      # Power coefficient - defined as .txt in "03_Turbine_Data"
    "Ct" => [3.0, 4.0],      # Thrust coefficient - defined as .txt in "03_Turbine_data"
    "Yaw" => [225, 225, 225], # Yaw angle of the turbines (In geographical DEG)

    ##########      (3) Atmospheric data       ######################
    #Use either 3.1 for single computation OR 3.2 for AEP computation
    # (3.1) Single computatiosn 
    #       This section is only used for single case computation    
    "u_ambient" => 10.0, # Ambient wind speed in [m/s]
    "alpha" => 225,    # Geographical direction of the wind speed in [°]. -> N == 0°
    "TI_a" => 0.07,    # Ambient turbulence intensity in [%]
    # (3.2) AEP computation 
    #       This section is only used for AEP computation      
    "Wind_rose" => 3.14159999, # Get wind rose as specified in "04_Ambient_data"

    ##########      (4) Computational setting  ######################
    "SimpleComp" => true,
    "AEPComp" => false,
    "Optim" => false,

    ##########      (5) Numerical parameters   ######################
    "RotorRes" => 100, #Number of points used to distrectisize the turbine's rotors

    ##########      (6) Graphical output       ######################
    "z" => 100

)
