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
        0.037617442870482, 0.70386928755986, 2.06862496038165, 2.96296058633423, 
        4.19112689408088, 4.9610745700036, 5.36590909310877, 5.99336081493984, 
        12.7388667463041, 8.04735461112754, 4.65874030506133, 12.5278425047245, 
        7.53213428352991, 12.9832434749864, 11.4439656097754, 15.6245817527353, 
        15.3221190536884, 21.6501206340332, 9.01786914438853, 14.2854422862743, 
        22.5733143394353, 19.9001279829305, 17.4201794757227, 9.53002429871616, 
        31.0093936523177, 20.1815101018876, 22.7061131702835, 17.7800056604986, 
        15.3757687434198, 29.0447584281646, 26.423965230646, 29.2407895303151, 
        36.5794923927964, 25.7738596012576, 27.1646353726186, 28.3159525083547, 
        31.5097305913561, 35.9116862033795, 32.4809647078191, 35.0108435099233, 
        40.4906399183072, 47.4910608013055, 56.7682258802251, 38.2379881777972, 
        32.4895845817412, 45.8704486847747, 39.9699889063481, 48.3601808109306, 
        46.5193813254214, 51.0675852336404, 41.7773717893322, 52.2270493564652, 
        34.7849354199214, 41.5078560950976, 47.7660073328075, 41.2370396214249, 
        52.8, 55.9816241943712, 60.2132787855731, 66.1543490075053, 
        51.6237053143583, 52.6451344485819, 61.9733145696716, 59.5323196036149, 
        58.7, 63.8543884597293, 64.5614315823471, 57.9172372361691, 
        56.6975992092866, 66.927416400415, 68.3161105079751, 69, 
        63.1501306083627, 37.8976029730991, 48.5969640909351, 65.3121955177782, 
        61.2131669543962, 67.678027356928, 54.0473078087759, 64.834553125265
    ], 
    "y_vec" => -[
        0.013484557804708, 4.73351577923936, 16.7584627674223, 23.5561910196013, 
        33.968147758905, 37.9950801243123, 43.4382325738325, 47.5427129721352, 
        1.52203982701976, 3.68169406758461, 12.4248564239154, 11.7375794487288, 
        25.3631837550851, 32.9052322419763, 40.5369829442608, 47.8837210046803, 
        7.18025042313192, 9.26341133970766, 21.9735849852595, 21.0384956784569, 
        28.4370965195701, 35.846209332528, 42.705757175383, 48.6299164297325, 
        1.99145905632608, 0.563710581584548, 17.1938370533637, 15.3498809920574, 
        26.7409286474901, 33.5495804504402, 40.6082315799248, 47.0578400397897, 
        1.14532648531499, 2.83888291727143, 10.8800681361102, 19.4178433623723, 
        27.3942531526216, 30.3174597893915, 42.3791781849196, 48.6289472561082, 
        0.640617936403942, 6.75977679761558, 11.349853974459, 14.810535195176, 
        12.6092245354078, 25.9119916445893, 41.7270646528813, 47.9656883336025, 
        0.228682899450648, 9.1574663810562, 13.6858375301102, 17.5229508183637, 
        21.8724182139961, 32.0799846983618, 41.155993074756, 45.5553097173971, 
        5.24e-18, 6.34389230910802, 14.1844046521093, 25.5659987217377, 
        27.9672987724243, 35.2646976922068, 42.0509551090612, 48.6299984701831, 
        2.08e-18, 8.26549497937587, 12.6554029615115, 19.6571385811674, 
        29.7516518417668, 33.1984299439432, 43.0870756842053, 48.63, 
        3.1706095256668, 4.83247295914046, 19.3937652045311, 18.7403447046287, 
        21.3384617608537, 38.7281382851855, 44.1013156639835, 47.8047852393555
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
