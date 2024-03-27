#= 
Script to Excecute & Initialise 
=#
module Preprocessing
export generateWF
include("PKG_Manager.jl") #make neccecary packages available (using XX)

# Function to generate struct instances from files in Input folder
function generateWF(prefix::AbstractString, folder::AbstractString)
    # Get a list of all files in the Inputfolder
    files = readdir(folder)
    # Filter out files that start with the given prefix
    script_files = filter(file -> startswith(file, prefix), files)
    # Find the number of applicable inputfiles
    n=length(script_files);   
    WF=Vector{Windfarm}(undef, n);  # Creating array to contain all input data
    
    # Execute each script -> Create array of mutable structs filled with information from the input file
    for i=1:n
        # User Input
        include(joinpath(folder, script_files[i]))
        # Create an instance of the struct
        WF[i] = WFConstructor(
            userdata["name"],
            userdata["N"],
            userdata["x_vec"],
            userdata["y_vec"],
            userdata["D"],
            userdata["H"],
            userdata["Cp"],
            userdata["Ct"],
            userdata["u_ambient"],
            userdata["alpha"],
            userdata["TI_a"],
            userdata["Wind_rose"],
            userdata["x"],
            userdata["y"],
            userdata["z"]
        )
    end
return WF
end

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
        x::Int;
    ##########      (5) Numerical parameters   ######################
        y::Int;
    ##########      (6) Graphical output       ######################
        z::Int;
end
    
# Constructor function to create instances of the struct Windfarm
function WFConstructor(name::String, N::Int,  x_vec::Vector{Float64}, y_vec::Vector{Float64}, D::Float64, H::Float64, Cp::Vector{Float64}, Ct::Vector{Float64}, u_ambient::Float64, alpha::Float64, TI_a::Float64, Wind_rose::Float64, x::Int, y::Int, z::Int)
        windfarm_instance = Windfarm(name, N,  x_vec, y_vec, D, H, Cp, Ct, u_ambient, alpha, TI_a, Wind_rose, x, y, z)
        return windfarm_instance
end

end