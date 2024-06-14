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
    "y_vec" => [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49],     #zeros(80,),         #Y-Coordinates [0, 0, 0, 0, 0, 0],

    ##########      (2) Turbine data           ######################
    "Yaw" => 90 .+ zeros(80,),     # Yaw angle of the turbines (In geographical DEG)
    "Turbine_Type" => "VestasV80",  #Turbine Type. One type for the whole wind farm. 
                                    #Possible Inputs: "VestasV80", "NREL_5MW", "DTU_10MW", "IEA_15MW"

    ##########      (3) Atmospheric data       ######################
    #Use either 3.1 for single computation OR 3.2 for AEP computation
    # (3.1) Single computatiosn 
    #       This section is only used for single case computation    
    "u_ambient" => 8,      # [m/s] Ambient wind speed
    "alpha"     => 90,     # [°] Geographical direction of the wind. -> N == 0°
    "TI_a"      => 0.045,    # [-] Ambient turbulence intensity in [-]
    "z_Surf"    => 0.005,   # [-] Surface roughness of the modelled case *for offshore conditions z_Surf should equal between 0.0001 (calm see) and 0.01 (high waves)
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
    #Correction Models
    "Meandering"=> true,                #Meandering correction as proposed by Braunbehrens & Segalini (2019).

    ##########      (5) Numerical parameters   ######################
    "Dimensions"            => "3D",        #Choose dimensions resolution. 1) Three dimensional space or 2) two dimansional plane at Hub height.
                                            #Possible inputs: "3D", "2D"
    "Rotor_Discretization"  => "gridded",   #Specifies the rotor descritization technique. Current choices: 1) Evenly distributed grid (slow with small error), 2) Fibonacci-Latice distributed points (quicker). 
                                            #Possible inputs: "gridded", "fibonacci" !!! Gridded has to be checked and corrected/ Thrown out
    "Rotor_Res"             => 100,          #Number of points used to represent the rotor. Reccomendation: 100 for "griddeed" & XX for "fibonacci".

    ##########      (6) Graphical output       ######################
    # Simple plots, no further computation:
    "Plot_power"        => true,    #Plots power output of several turbines
    "Plot_windspeed"    => true,   #Plots average inflow windspeed of several turbines
    "Plot_turbulence"   => false,   #Plots average inflow turbulence of several turbines
    "Turbine_Identification"    => [41, 42, 43, 44, 45, 46, 47, 48, 49, 50], #Identify, which turbines should be included in the plot
    "Normalize_to"              => 41,                               #Specify which turbines power the plot should be normalised to (If no normalisation is wanted, type: 0)
    # Advanced plots, advanced computation will commence
    "Plot_wind_field"       => false,   #Plots wind field for one simple case
    "Plot_turbulence_field" => false,   #Plots turbulence field for one simple case
    "Wind_Direction"        => 270.0,     #Wind direction for plot (has to be a direction included during computation!)


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
