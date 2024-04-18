#= 
Module for:
    1) Reading input data
    2) Initialising arrays & allocate space for computation
=#
module Input_Processing
using Pkg, FileIO, Parameters, InteractiveUtils, DataStructures, Revise

export generateWF, ComputationData

# Function to generate struct instances from files in Input folder
function generateWF(prefix::AbstractString, folder::AbstractString)
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
        WF[i] = WFConstructor(userdata)
    end
return WF
end #generateWF
    
# Constructor function to create instances of the struct Windfarm
function WFConstructor(userdata::OrderedDict{String, Any})
    args = (userdata[arg] for arg in keys(userdata)) #Get all fields within the dictionary
    return Windfarm(args...)                         #Pass them one by one to the strtuct definition
end #WFConstructor

# Initiate template struct, generation function.
# Initialises all variables to be assigned from the input files.
mutable struct Windfarm
    ##########      (1) Wind farm data         ######################
        # Name of the wind Farm
        name::String;
        
        # Wind Farm Data
        N::Int; #Number of turbines
        
        x_vec::Vector{Float64};   #X-Coordinates
        y_vec::Vector{Float64};        #Y-Coordinates
    
    ##########      (2) Turbine data           ######################
    
        D::Float64;             # Turbine diameter in [m]
        H::Float64;             # Hub height in [m]
        Cp::Vector{Float64};    # Power coefficient - defined as .txt in "03_Turbine_Data"
        Ct::Vector{Float64};    # Thrust coefficient - defined as .txt in "03_Turbine_data"
        Yaw::Vector{Float64};   # Yaw angle of each turbine
    
    ##########      (3) Atmospheric data       ######################
    #Use either 3.1 for single computation OR 3.2 for AEP computation
    
    # (3.1) Single computation 
    #       This section is only used for single case computation    
        u_ambient::Float64; # Ambient wind speed in [m/s]
        alpha::Float64;     # Geographical direction of the wind speed in [°]. -> N == 0°
        TI_a::Float64;      # Ambient turbulence intensity in [%]
    
    # (3.2) AEP computation 
    #       This section is only used for AEP computation  #   
        
        Wind_rose::Float64; # Get wind rose as specified in "04_Ambient_data"
    
    ##########      (4) Computational setting  ######################
        SimpleComp::Bool;
        AEPComp::Bool;
        Optimisation::Bool;
    ##########      (5) Numerical parameters   ######################
        RotorRes::Int;
    ##########      (6) Graphical output       ######################
        z::Int;
end # mutable struct Windfarm

end #module