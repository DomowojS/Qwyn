# Creates a dictionary of Input data
using JLD2 #Needed for loading turbine & meteorological data
using LinearAlgebra
userdata = OrderedDict{String, Any}(
    ##########      (1) Wind farm data         ######################
    # Name of the wind Farm
    "name" => "WF_Evaluate_Optim",
    # Wind Farm Data
    "N" => 80, #Number of turbines  
    "x_vec" => [0.06738891, 1.099614562, 2.02675399, 3.967611421, 2.997937965, 4.091642564, 4.8458306, 
    5.916883429, 11.05438351, 7.026528672, 9.136814506, 4.580068629, 7.167099305, 9.847888817, 
    11.35745694, 7.169343814, 14.72055849, 15.80730831, 17.23833311, 20.37407164, 16.7199139, 
    13.93482122, 21.47104727, 18.0768328, 27.09091375, 18.77225499, 19.03225921, 22.52556038, 
    21.8135041, 27.53019447, 28.8995724, 25.03739024, 40.07591255, 25.19924888, 28.66032984, 
    31.09074226, 34.06851861, 35.28949372, 36.49358861, 30.67907845, 53.95047862, 37.71096089, 
    27.19506942, 40.87212139, 34.61186175, 48.91241393, 42.31540327, 37.53700157, 47.48645533, 
    60.47131697, 32.54334507, 59.76704873, 46.02952669, 63.1508116, 46.13691112, 53.83341752, 
    64.62072854, 55.29823942, 44.61058299, 49.79374158, 51.43099818, 54.61712543, 59.36972586, 
    69, 63.83244125, 43.4420702, 62.04423174, 65.47359114, 66.8964066, 57.52131586, 61.21192394, 
    64.73727513, 62.87585351, 49.04621458, 56.98907474, 38.44769745, 66.14483603, 67.46678565, 
    52.01825794, 68.18127126
    ], 
    "y_vec" => -[0, 6.069910461, 13.10593258, 18.18164994, 24.25800656, 33.13245139, 39.24824042, 
    44.27367739, 1.16E-17, 2.064903688, 16.12582858, 7.950094185, 26.15182612, 36.06108148, 
    42.57164411, 48.63, 2.865115486, 6.972728126, 12.86471982, 23.54203342, 31.22351293, 
    37.93539483, 40.8685401, 48.63, 7.14E-17, 1.723950604, 18.8638856, 10.70486735, 
    28.3554422, 30.36124241, 43.65789193, 48.62816713, 0.117247823, 4.457145668, 14.76647897, 
    20.94986449, 29.27964137, 34.6602265, 40.51995019, 48.61324797, 4.09E-17, 4.802914801, 
    9.257809774, 14.29080387, 11.81580762, 33.93783581, 44.31377229, 47.06185195, 0.008572025, 
    4.574580102, 1.572221221, 22.82513236, 25.36323285, 33.22049679, 42.55787221, 48.62980252, 
    13.13600479, 5.574694999, 16.08746372, 18.46056958, 27.25049126, 37.49875786, 45.82723386, 
    48.63, 6.746936301, 7.229515292, 17.86959666, 20.04845616, 31.58037552, 29.46941651, 
    39.81562298, 47.4671058, 0.004442421, 9.075599699, 14.47259316, 23.2181557, 25.48889602, 
    36.20329772, 42.20016632, 41.99420359
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
    "Meandering"=> true,                #Meandering correction as proposed by Braunbehrens & Segalini (2019).

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
