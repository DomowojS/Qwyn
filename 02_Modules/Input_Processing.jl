#= 
Module for:
    1) Reading input data
    2) Initialising arrays & allocate space for computation
=#
module Input_Processing
using Pkg, FileIO, Parameters, InteractiveUtils, DataStructures, Revise
export generateWF, ComputationData

    function generateWF(prefix::AbstractString, folder::AbstractString, u_ambient::Real, alpha::Real, TI_a::Real)
    #Function to generate struct instances from files in Input folder
        # Get a list of all files in the Inputfolder
        files = readdir(folder)
        # Filter out files that start with the given prefix
        script_files = filter(file -> startswith(file, prefix), files)
        # Find the number of applicable inputfiles
        n=length(script_files);   
        WF=Vector{Windfarm}(undef, n);  # Creating array to contain all input data (of the type of Windfarm - mutable struct)

            # Iterate over all input files and store them in struct array
            for i=1:n
                # User Input
                include(joinpath(folder, script_files[i])) #Overwrite user data after each iteration
                # Create an instance of the struct
                WF[i] = WFConstructor(u_ambient, alpha, TI_a, userdata)
            end

        return WF
    end #generateWF
    
    function WFConstructor(u_ambient::Real, alpha::Real, TI_a::Real, userdata::OrderedDict{String, Any})
        #Placeholders
        userdata["D"]         = 0.0                     # Turbine diameter 
        userdata["H"]         = 0.0                     # Hub height in [m]
        userdata["P_Input"]   = zeros(1,1)              # Power coefficient - defined as .txt in "03_Turbine_Data"
        userdata["Ct_Input"]  = zeros(1,1)              # Thrust coefficient - defined as .txt in "03_Turbine_data"
        userdata["u_ambient_zprofile"] = zeros(1,1,1)   # [m/s] height profile of the wind as vector of z coordinates resulting from amount of rotor resolution points 
        #UserInput
        userdata["u_ambient"] = u_ambient               # [m/s] Ambient wind speed
        userdata["alpha"]     = alpha                   # [°] Geographical direction of the wind. -> N == 0°
        userdata["TI_a"]      = TI_a                    # [-] Ambient turbulence intensity in [-]
        # Extract all fields within the updated dictionary
        args = (userdata[arg] for arg in keys(userdata))

        # Return the Windfarm struct by passing all arguments
        return Windfarm(args...)
    end


    mutable struct Windfarm
    #Input struct definition
        ##########      (1) Wind farm data         ######################
            # Name of the wind Farm
            name::String;
            
            # Wind Farm Data
            N::Int; #Number of turbines
            
            x_vec::Vector{Float64};   #X-Coordinates
            y_vec::Vector{Float64};   #Y-Coordinates
        
        ##########      (2) Turbine data           ######################   
            Yaw::Vector{Float64};   # Yaw angle of each turbine
            Turbine_Type::String;        # Turbine Type

        ##########      (3) Atmospheric data       ######################
        #Use either 3.1 for single computation OR 3.2 for AEP computation
        
        # (3.1) Single computation 
        #       This section is only used for single case computation    
            z_Surf::Float64;    # [-] Surface roughness of the modelled case *for offshore conditions z_Surf should equal between 0.0001 (calm see) and 0.01 (high waves)
            z_r::Float64;       # [m] Height the average wind speed "u_ambient" was measured. If not known, choose z_u = 10
        
        # (3.2) AEP computation 
        #       This section is only used for AEP computation  #   
            
            Wind_rose::Float64; # Get wind rose as specified in "04_Ambient_data"
        
        ##########      (4) Computational setting  ######################
            SimpleComp::Bool;           # For the computation of one case
            AEPComp::Bool;              # For the estimation of the farms AEP
        ## Advanced settings:
        #Superposition Method
            Superpos::String;    # Superposition using linear rotorbased summation for velocity deficit
        #Correction Models
            Meandering::Bool;    #Meandering correction as proposed by Braunbehrens & Segalini (2019).

        ##########      (5) Numerical parameters   ######################
            Rotor_Discretization::String; #Specifies the rotor descritization technique. Current choices: 1) Evenly distributed grid (slow with small error), 2) Fibonacci-Latice distributed points (quicker). Possible inputs: "gridded", "fibonacci"
            Rotor_Res::Int;               #Number of points used to represent the rotor. Reccomendation: 100 for "griddeed" & XX for "fibonacci".
            # For Momentum conserving superposition only:
            Uc_Res::Int;                  #Number of points to comppute wake for convection velocity

        ##########      (6) Graphical output       ######################
            # Simple plots, no further computation:
            Plot_power::Bool;       #Plots power output of several turbines
            Plot_windspeed::Bool;   #Plots average inflow windspeed of several turbines
            Plot_turbulence::Bool;  #Plots average inflow turbulence of several turbines
            Turbine_Identification::Vector{Int};    #Identify, which turbines should be included in the plot
            Normalize_to::Int;                      #Specify which turbines power the plot should be normalised to (If no normalisation is wanted, type: 0)
            # Advanced plots, advanced computation will commence
            Plot_wind_field::Bool;      #Plots wind field for one simple case
            Plot_turbulence_fiel::Bool; #Plots turbulence field for one simple case
            Wind_Direction::Float64;    #Wind direction for plot (has to be a direction included during computation!)

            z::Int;    
            Z_Max::Float64;  #Maximum height
            Z_Min::Float64;  #Minimum height

        ##########   Turbine Placeholders           ######################
            D::Float64;             # Turbine diameter in [m]
            H::Float64;             # Hub height in [m]
            P_Input::Matrix{Float64};    # Power coefficient - defined as .txt in "03_Turbine_Data"
            Ct_Input::Matrix{Float64};    # Thrust coefficient - defined as .txt in "03_Turbine_data"
            #Atmospheric data placeholders:
            u_ambient_zprofile::Array{Float64,3}; # [m/s] height profile of the wind as vector of z coordinates resulting from amount of rotor resolution points 
        ##########   Direct User Input              ######################
            u_ambient::Float64; # [m/s] Ambient wind speed
            alpha::Real;     # [°] Geographical direction of the wind speed. -> N == 0°
            TI_a::Real;      # [-] Ambient turbulence intensity in [-]
    end # mutable struct Windfarm

end #module