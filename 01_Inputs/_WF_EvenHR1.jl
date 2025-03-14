using JLD2 #Needed for loading turbine & meteorological data
    using LinearAlgebra
    userdata = OrderedDict{String, Any}(
        ##########      (1) Wind farm data         ######################
        # Name of the wind Farm
        "name" => "Even_HR1",
        # Wind Farm Data
        "N" => 80, #Number of turbines  
        "x_vec" => [0,	7,	14,	21,	28,	35,	42,	49,	0,	7,	14,	21,	28,	35,	42,	49,	0,	7,	14,	21,	28,	35,	42,	49,	0,	7,	14,	21,	28,	35,	42,	49,	0,	7,	14,	21,	28,	35,	42,	49,	0,	7,	14,	21,	28,	35,	42,	49,	0,	7,	14,	21,	28,	35,	42,	49,	0,	7,	14,	21,	28,	35,	42,	49,	0,	7,	14,	21,	28,	35,	42,	49,	0,	7,	14,	21,	28,	35,	42,	49],
        "y_vec" => [0,	0,	0,	0,	0,	0,	0,	0,	7,	7,	7,	7,	7,	7,	7,	7,	14,	14,	14,	14,	14,	14,	14,	14,	21,	21,	21,	21,	21,	21,	21,	21,	28,	28,	28,	28,	28,	28,	28,	28,	35,	35,	35,	35,	35,	35,	35,	35,	42,	42,	42,	42,	42,	42,	42,	42,	49,	49,	49,	49,	49,	49,	49,	49,	56,	56,	56,	56,	56,	56,	56,	56,	63,	63,	63,	63,	63,	63,	63,	63],

        ##########      (2) Turbine data           ######################
        "Yaw" => 270 .+ zeros(80,),     # Yaw angle of the turbines (In geographical DEG)
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
        ## (4.1) Advanced Computational setting:
        #Wake Model
        "WakeModel" => "Ishihara", #Single wake model. Choose between Ishihara-Qian (2018) and TurbOPark (2022). 
                                   #Possible inputs: "Ishihara", "TurbOPark"
    
        #Superposition Method
        "Superpos"  => "Linear_Rotorbased", #Superposition method for velocity deficits. Choose between linear rotorbased summation & momentum conserving approach. 
                                              #Possible inputs: "Linear_Rotorbased", "Momentum_Conserving"
        #Correction Models
        "Meandering"=> false,                #Meandering correction as proposed by Braunbehrens & Segalini (2019).
    
        ##########      (5) Numerical parameters   ######################
        "Rotor_Discretization"  => "smart_grid",#Specifies the rotor descritization technique. Current choices: 1) Evenly distributed grid (slow with small error), 2) Fibonacci-Latice distributed points (quicker). 3) Smart grid. Mixture of Gaussian qudrature/ Circular gauss and fibonacci
                                                #Possible inputs: "gridded", "fibonacci", "smart_grid" !!! Gridded has to be checked and corrected/ Thrown out
        "Rotor_Res"             => 12,           #Number of points used to represent the rotor. Reccomendation: 100 for "gridded" & >21 for "fibonacci".
                                                #For Smart Grid, Possible Inputes are: "1, 4, 6, 7, 9, 12, 21, >21".
            # For Momentum conserving superposition only:
            "Uc_Res"                => 1000,      #Number of points to comppute wake for convection velocity (at each relevant streamwise position x). Has to be > 4
    
        ##########      (6) Result struct request  ######################
        
        "Extended_Output"   => true,    #"false" returns consice result struct with the most important input and computed turbine performance data.
                                        #"true" returns all input computation struct in cell array "WF" as well as full "Computation_Struct" which includes all computational arrays & results. 
        
        ##########      (7) Graphical output       ######################
        # Simple plots, no further computation:
        "Plot_power"        => false,    #Plots power output of several turbines
        "Plot_windspeed"    => false,  #Plots average inflow windspeed of several turbines
        "Plot_turbulence"   => false,   #Plots average inflow turbulence of several turbines
        "Turbine_Identification"    => [1, 2, 3], #Identify, which turbines should be included in the plot
        "Normalize_to"              => 1,         #Specify which turbines power the plot should be normalised to (If no normalisation is wanted, type: 0)
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
    
