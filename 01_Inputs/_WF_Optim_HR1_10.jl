# Creates a dictionary of Input data
using JLD2 #Needed for loading turbine & meteorological data
using LinearAlgebra
userdata = OrderedDict{String, Any}(
    ##########      (1) Wind farm data         ######################
    # Name of the wind Farm
    "name" => "WF_Evaluate_Optim",
    # Wind Farm Data
    "N" => 80, #Number of turbines  
    "x_vec" => [
        0.126478, 0.886828526, 1.606310513, 3.151883896, 3.889586006, 
        4.592078142, 5.32969164, 6.224295212, 14.71239207, 6.326532793, 
        8.58892905, 7.265318323, 9.94029514, 15.01084363, 9.941838943, 
        16.50801979, 4.305634633, 11.3678837, 12.42255688, 13.54075812, 
        20.91578897, 19.25062424, 15.70544238, 21.7543806, 34.1011339, 
        25.6089909, 22.58151478, 23.42085377, 29.41034497, 30.5651173, 
        25.80472061, 26.86845776, 39.67912178, 15.78054908, 29.82562716, 
        32.89618506, 35.01189578, 44.8737627, 40.41635619, 38.79438384, 
        48.95534314, 28.31363152, 35.72582515, 52.28655445, 22.06163269, 
        42.56059226, 46.12541955, 50.4723078, 57.47025606, 41.89429964, 
        38.72838976, 45.54466096, 54.41876941, 38.59376669, 49.79347618, 
        56.83405553, 51.03162876, 59.90897017, 64.96664782, 44.7171219, 
        56.32044418, 59.22545066, 60.08788787, 62.9440941, 62.9213327, 
        64.25147339, 61.14220051, 65.58972176, 66.81518703, 63.10626097, 
        64.91054445, 33.88295283, 63.62685614, 55.58474454, 37.60777165, 
        48.89159216, 66.16550691, 67.44908225, 67.93251306, 68.83622216
    ], 
    "y_vec" => -[
        0.693487619, 6.963033381, 12.24063574, 25.54305451, 29.32078272, 
        33.60144143, 41.38823043, 48.63, 0, 5.332027675, 23.1418028, 
        16.74187491, 31.14671234, 35.11775993, 43.67056029, 48.63, 
        9.51E-05, 5.023641354, 14.49288862, 20.96790233, 27.68728298, 
        32.85059678, 43.07684448, 48.03172764, 0.580070541, 1.555686741, 
        16.54892257, 22.21325285, 26.71886386, 36.5668688, 42.41606083, 
        47.33476987, 0.000158554, 6.814058487, 13.02885163, 20.76946633, 
        25.93161066, 34.97179307, 43.59369609, 48.63, 3.968121898, 
        8.538097865, 5.312977028, 10.83134142, 6.382537892, 28.44536046, 
        39.60715109, 48.63, 0.14916021, 6.784502646, 16.64796539, 
        18.95771683, 22.28425074, 33.51530362, 42.71141453, 48.30296598, 
        0.000261846, 6.648845327, 16.00065958, 14.09950569, 26.97499435, 
        34.11406933, 40.84848623, 47.90451319, 0, 10.14319179, 
        17.72494603, 20.98969484, 30.95421446, 32.58070303, 44.50364697, 
        46.41422247, 5.080669001, 3.368192153, 12.01432079, 20.53507273, 
        25.65643351, 36.05981167, 39.97801839, 47.30258065
    ],

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
    "Rotor_Res"             => 10,           #Number of points used to represent the rotor. Reccomendation: 100 for "gridded" & XX for "fibonacci".
    
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
    "Normalize_to"              => 4,                                       #Specify which turbines power the plot should be normalised to (If no normalisation is wanted, type: 0)
    # 2D surface plots, advanced computation will commence
    "Plot_wind_field"       => false,   #Plots wind field for one simple case
    "Plot_turbulence_field" => false,   #Plots turbulence field for one simple case
    "Wind_Direction"        => 270.0,   #Wind direction for plot (has to be a direction included during computation!)
    "Resolution"            => 0.5,     #Pick resolution in terms of Diameter  
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
