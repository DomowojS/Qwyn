# Creates a dictionary of Input data
using JLD2 #Needed for loading turbine & meteorological data
using LinearAlgebra
userdata = OrderedDict{String, Any}(
    ##########      (1) Wind farm data         ######################
    # Name of the wind Farm
    "name" => "3Turb",
    # Wind Farm Data
    "N" => 3, #Number of turbines  
    "x_vec" => [0, 4, 18],#[0, 0.857000000000000, 1.71400000000000, 2.57000000000000, 3.43000000000000, 4.29000000000000, 5.14000000000000, 6, 7, 7.85700000000000, 8.71400000000000, 9.57000000000000, 10.4300000000000, 11.2900000000000, 12.1400000000000, 13, 14, 14.8570000000000, 15.7140000000000, 16.5700000000000, 17.4300000000000, 18.2900000000000, 19.1400000000000, 20, 21, 21.8570000000000, 22.7140000000000, 23.5700000000000, 24.4300000000000, 25.2900000000000, 26.1400000000000, 27, 28, 28.8570000000000, 29.7140000000000, 30.5700000000000, 31.4300000000000, 32.2900000000000, 33.1400000000000, 34, 35, 35.8570000000000, 36.7140000000000, 37.5700000000000, 38.4300000000000, 39.2900000000000, 40.1400000000000, 41, 42, 42.8570000000000, 43.7140000000000, 44.5700000000000, 45.4300000000000, 46.2900000000000, 47.1400000000000, 48, 49, 49.8570000000000, 50.7140000000000, 51.5700000000000, 52.4300000000000, 53.2900000000000, 54.1400000000000, 55, 56, 56.8570000000000, 57.7140000000000, 58.5700000000000, 59.4300000000000, 60.2900000000000, 61.1400000000000, 62, 63, 63.8570000000000, 64.7140000000000, 65.5700000000000, 66.4300000000000, 67.2900000000000, 68.1400000000000, 69], #[0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 0, 7, 14, 21, 28, 35, 42, 49, 56, 63], #LinRange(0,21,80),   #X-Coordinates [0, 7, 14, 21, 28, 35],
    "y_vec" => [0, 0, 0],#-[0, 6.95000000000000, 13.8900000000000, 20.8400000000000, 27.7890000000000, 34.7370000000000, 41.6800000000000, 48.6300000000000, 0, 6.95000000000000, 13.8900000000000, 20.8400000000000, 27.7890000000000, 34.7370000000000, 41.6800000000000, 48.6300000000000, 0, 6.95000000000000, 13.8900000000000, 20.8400000000000, 27.7890000000000, 34.7370000000000, 41.6800000000000, 48.6300000000000, 0, 6.95000000000000, 13.8900000000000, 20.8400000000000, 27.7890000000000, 34.7370000000000, 41.6800000000000, 48.6300000000000, 0, 6.95000000000000, 13.8900000000000, 20.8400000000000, 27.7890000000000, 34.7370000000000, 41.6800000000000, 48.6300000000000, 0, 6.95000000000000, 13.8900000000000, 20.8400000000000, 27.7890000000000, 34.7370000000000, 41.6800000000000, 48.6300000000000, 0, 6.95000000000000, 13.8900000000000, 20.8400000000000, 27.7890000000000, 34.7370000000000, 41.6800000000000, 48.6300000000000, 0, 6.95000000000000, 13.8900000000000, 20.8400000000000, 27.7890000000000, 34.7370000000000, 41.6800000000000, 48.6300000000000, 0, 6.95000000000000, 13.8900000000000, 20.8400000000000, 27.7890000000000, 34.7370000000000, 41.6800000000000, 48.6300000000000, 0, 6.95000000000000, 13.8900000000000, 20.8400000000000, 27.7890000000000, 34.7370000000000, 41.6800000000000, 48.6300000000000],

    ##########      (2) Turbine data           ######################
    "Yaw" => 270 .+ zeros(3,),     # Yaw angle of the turbines (In geographical DEG)
    "Turbine_Type" => "VestasV80",  #Turbine Type. One type for the whole wind farm. 
                                    #Possible Inputs: "VestasV80", "NREL_5MW", "DTU_10MW", "IEA_15MW"

    ##########      (3) Ambient data       ######################
    #Use either 3.1 for single computation OR 3.2 for AEP computation
    # (3.1) Single computatiosn 
    #       This section is only used for single case computation    
    "z_Surf"    => 0.001,   # [-] Surface roughness of the modelled case *for offshore conditions z_Surf should equal between 0.0001 (calm see) and 0.01 (high waves)
    "z_r"       => 70.0,    # [m] Height the average wind speed "u_ambient" was measured. If not known, choose z_r = 10

    # (3.2) AEP computation 
    #       This section is only used for AEP computation      
    "Wind_rose" => 3.14159999, # Get wind rose as specified in "04_Ambient_data"

    ##########      (4) Computational setting  ######################
    "CompSetting"    => "Simple",    #= What do you want to compute? 1: "Simple" - comp. of one case
                                                                     2: "AEP"    - comp of AEP  
                                                                     3: "Optim"  - WF Optimization =#
    
    ## (4.1) Advanced Computational setting:
    #Superposition Method
    "Superpos"  => "Linear_Rotorbased", #Superposition method for velocity deficits. Choose between linear rotorbased summation & momentum conserving approach. 
                                        #Possible inputs: "Linear_Rotorbased", "Momentum_Conserving"
    #Correction Models
    "Meandering"=> false,                #Meandering correction as proposed by Braunbehrens & Segalini (2019).

    ##########      (5) Numerical parameters   ######################
    "Rotor_Discretization"  => "fibonacci", #Specifies the rotor descritization technique. Current choices: 1) Evenly distributed grid (slow with small error), 2) Fibonacci-Latice distributed points (quicker). 
                                            #Possible inputs: "gridded", "fibonacci" !!! Gridded has to be checked and corrected/ Thrown out
    "Rotor_Res"             => 4,           #Number of points used to represent the rotor. Reccomendation: 100 for "gridded" & XX for "fibonacci".
    
        # For Momentum conserving superposition only:
        "Uc_Res"                => 1000,      #Number of points to comppute wake for convection velocity (at each relevant streamwise position x). Has to be > 4

    ##########      (6) Result struct request  ######################
    
    "Extended_Output"   => false,    #"false" returns consice result struct with the most important input and computed turbine performance data.
                                    #"true" returns all input computation struct in cell array "WF" as well as full "Computation_Struct" which includes all computational arrays & results. 
    
    ##########      (7) Graphical output       ######################
    # Simple plots, no further computation:
    "Plot_power"        => false,    #Plots power output of several turbines
    "Plot_windspeed"    => false,   #Plots average inflow windspeed of several turbines
    "Plot_turbulence"   => false,   #Plots average inflow turbulence of several turbines
    "Turbine_Identification"    => [4, 12, 20, 28, 36, 44, 52, 60, 68, 76], #Identify, which turbines should be included in the plot
    "Normalize_to"              => 4,                                      #Specify which turbines power the plot should be normalised to (If no normalisation is wanted, type: 0)
    # 2D surface plots, advanced computation will commence
    "Plot_wind_field"       => true,   #Plots wind field for one simple case
    "Plot_turbulence_field" => false,   #Plots turbulence field for one simple case
    "Wind_Direction"        => 270.0,   #Wind direction for plot (has to be a direction included during computation!)
    "Resolution"            => 0.1,     #Pick resolution in terms of Diameter  
    "Height"                => 70.0,    #At what height [m] do you want to plot the 2D wind field (Y-Z plane)?
    "Depth"                 => 0.0,     #At what Y coordinate [1/D] do you want to plot the cross section wind field (X-Z plane)?
    ## Full 3D plot lots of RAM required. Uses same resolution as 2D surface plots
    #"Plot_wind_field_3D"       => false,   #Plots wind field for one simple case
    #"Plot_turbulence_field_3D" => false,   #Plots turbulence field for one simple case

    ##########   Literature Input              ######################
    #= These Numbers are placeholders. They overwritten by data from literature.
    You can provide this as .jld2 file for turbine & atmospheric data.
        1) Turbine data     --> provide .jld2 file providing turbine data in folder "04_Turbine_Data"
        2) Atmospheric data --> provide .jld2 file providing atmospheric data in folder "05_Atmospheric_Data"
    A documentation on how to provide the data can be found in the corresponding folders.
    =#
   )
