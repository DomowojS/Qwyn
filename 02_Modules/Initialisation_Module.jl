# Module to initialise matrices and compute all necessary coordinates

module Initialisation_Module
using JLD2, Interpolations, LinearAlgebra#, MAT, LatinHypercubeSampling, PlotlyJS, Colors

export initCompArrays, LoadTurbineDATA!, LoadAtmosphericData!, generate_rotor_grid

function initCompArrays(WindFarm)
    # Adjust user input for computation
    alpha_Comp  =   deg2rad(270 - WindFarm.alpha) #User Input logic: Geographical. Computation logic: Flow from left to right (270°(Western wind)==0°)
    Yaw_Comp    =   deg2rad.(270 .- WindFarm.Yaw)

    # Transform turbine positions according to inflow direction (Computation logic: windflow always left to right)
    TMPX = WindFarm.x_vec
    WindFarm.x_vec =  TMPX .* cos(alpha_Comp) .+ WindFarm.y_vec .* sin(alpha_Comp);
    WindFarm.y_vec = -TMPX .* sin(alpha_Comp) .+ WindFarm.y_vec .* cos(alpha_Comp);
    
    # Generates grid to represent the rotor plane
    if WindFarm.Rotor_Discretization == "gridded"
        #Generates an evenly distributed grid. Accurate representation for high number of points but very slow.
        Y_vec, Z_vec = generate_rotor_grid(WindFarm.Rotor_Res)
    elseif WindFarm.Rotor_Discretization == "fibonacci"
        #Generates a fibonacci-lattice distributed grid. Accurate representation for significantly lower amount of points.
        Y_vec, Z_vec = generate_fibonacci(WindFarm.RotorRes)
    else
        error("ERROR: Wrong choice of rotor discretization function for", WindFarm.Name,". Check variable Rotor_Discretization. Currently allowed entries: gridded, fibonacci.")
    end
    # Find and stor actual amount of points per rotor (after grid generation)
    Real_Rotor_Res = length(Y_vec);

    ### Compute coordinates (in D) for wake computation ###
    # Dimenisons: Relative X Coordinate, Relative Y Coordinate, Z Coordinate, Absolute turbine number
    XCoordinate = zeros(Float64, WindFarm.N , 1, WindFarm.N);  #Array for X coordinates of all points
    YCoordinate = zeros(Float64, WindFarm.N , Real_Rotor_Res, WindFarm.N);  #Array for Y Coordinates of all points
    ZCoordinate = zeros(Float64, 1 , Real_Rotor_Res, 1);  #Vector containing all height coordinates

    for i in 1:WindFarm.N
    # Create coordinate array of the structure: 1.Dim: Relative turbine, 2.Dim: RotorPoints, 3.Dim: Absolute turbine
        XCoordinate[:, :, i] .= XCoordinate[:, :, i] .+ WindFarm.x_vec .- WindFarm.x_vec[i];
        for j in 1:WindFarm.N
            YCoordinate[j, :, i] .= Y_vec .+ WindFarm.y_vec[j];
        end
        YCoordinate[:,:,i] .= YCoordinate[:,:,i] .- WindFarm.y_vec[i]
    end
    ZCoordinate[1,:,1] .= Z_vec

    #= Transform with respect to yaw angle
    TBDone!!
    =#


    # Create struct which holds all computation arrays
    CS=ComputationStruct(XCoordinate, YCoordinate, ZCoordinate, zeros(WindFarm.N , Real_Rotor_Res, WindFarm.N), Real_Rotor_Res,  alpha_Comp, Yaw_Comp,
                            zeros(1,1,WindFarm.N), zeros(1,1,WindFarm.N), 0, 0, (zeros(1,1,WindFarm.N) .+ WindFarm.u_ambient), zeros(1,1,WindFarm.N), (zeros(1,1,WindFarm.N) .+ WindFarm.TI_a),
                            zeros(1,1,WindFarm.N), zeros(1,1,WindFarm.N), zeros(1,1,WindFarm.N), zeros(1,1,WindFarm.N), zeros(1,1,WindFarm.N), 
                            zeros(1,1,WindFarm.N), zeros(1,1,WindFarm.N), zeros(1,1,WindFarm.N), zeros(WindFarm.N,Real_Rotor_Res,WindFarm.N), zeros(WindFarm.N,Real_Rotor_Res,WindFarm.N), 
                            zeros(WindFarm.N,Real_Rotor_Res,WindFarm.N), zeros(WindFarm.N,Real_Rotor_Res,WindFarm.N), zeros(1,Real_Rotor_Res,1), zeros(WindFarm.N,Real_Rotor_Res,WindFarm.N),
                            similar(XCoordinate, Bool), zeros(WindFarm.N,Real_Rotor_Res,WindFarm.N), zeros(WindFarm.N,Real_Rotor_Res,WindFarm.N), 100
                        )
 
    return WindFarm, CS

end #initCompArrays

function generate_rotor_grid(totalPoints::Int)
# This function generates a grid of points int the Y-Z plane to represent the turbine's rotor's.

    # Calculate approximate number of points per axis
    numPointsPerAxis = ceil(Int, sqrt(totalPoints / π))
    
    # Create a sufficiently dense grid to ensure we get the desired number of points
    x = range(-0.5, stop=0.5, length=numPointsPerAxis * 2)
    y = range(-0.5, stop=0.5, length=numPointsPerAxis * 2)
    X, Y = [xi for xi in x, yi in y], [yi for xi in x, yi in y]

    # Flatten the grid matrices
    X = vec(X)
    Y = vec(Y)

    # Calculate distance from the origin to each point
    distances = sqrt.(X.^2 .+ Y.^2)

    # Filter out points that are outside the rotorarea of radius 0.5
    insiderotor = distances .<= 0.5

    # Keep only the points inside the circle
    Y = Y[insiderotor]
    Z = X[insiderotor]

    return Y, Z
end

function LoadTurbineDATA!(WindFarm, CS)
#= This function loads all turbine data necessary for the computation
 It returns an update WindFarm & CS struct with Thrust coefficient and power curve, turbine diameter & hubheight. 
The ZCoordinate is also coorrected to have its origin at the Hubheigt of the turbine chosen for modelling. =#
    if WindFarm.Turbine_Type=="VestasV80"
        WindFarm.D = 80;
        WindFarm.H = 70;
        data = load("04_Turbine_Data\\VestasV80_2MW.jld2")  # Load data from Input 
        WindFarm.Ct_Input = data["CT_Input"];               # Thrustcoefficient vs. Windspeed
        WindFarm.P_Input  = data["P_Input"];                 # Power vs. Windspeed
        CS.Interp_Ct   =   LinearInterpolation(WindFarm.Ct_Input[:,1], WindFarm.Ct_Input[:,2]) #Interpolation function for Ct
        CS.Interp_P    =   LinearInterpolation(WindFarm.P_Input[:,1], WindFarm.P_Input[:,2]) #Interpolation function for P
        CS.Ct_vec   .=  CS.Interp_Ct(WindFarm.u_ambient);  #Ct of each turbine
        CS.P_vec    .=  CS.Interp_P(WindFarm.u_ambient);   #P of each turbine 
    elseif WindFarm.Turbine_Type=="NREL_5MW"
        WindFarm.D = 126;
        WindFarm.H = 90;
        data=load("04_Turbine_Data\\NREL_5MW.jld2")         # Load data from Input 
        WindFarm.Ct_Input;                                  # Thrustcoefficient vs. Windspeed
        WindFarm.P_Input;                                   # Power vs. Windspeed
        CS.Interp_Ct   =   LinearInterpolation(WindFarm.Ct_Input[:,1], WindFarm.Ct_Input[:,2]) #Interpolation function for Ct
        CS.Interp_P    =   LinearInterpolation(WindFarm.P_Input[:,1], WindFarm.P_Input[:,2]) #Interpolation function for P
        CS.Ct_vec   .=  CS.Interp_Ct(WindFarm.u_ambient);  #Ct of each turbine
        CS.P_vec    .=  CS.Interp_P(WindFarm.u_ambient);   #P of each turbine
    elseif WindFarm.Turbine_Type=="DTU_10MW"
        WindFarm.D = 119;
        WindFarm.H = 178.3;
        data=load("04_Turbine_Data\\DTU_10MW.jld2")         # Load data from Input 
        WindFarm.Ct_Input;                                  # Thrustcoefficient vs. Windspeed
        WindFarm.P_Input;                                   # Power vs. Windspeed
        CS.Interp_Ct   =   LinearInterpolation(WindFarm.Ct_Input[:,1], WindFarm.Ct_Input[:,2]) #Interpolation function for Ct
        CS.Interp_P    =   LinearInterpolation(WindFarm.P_Input[:,1], WindFarm.P_Input[:,2]) #Interpolation function for P
        CS.Ct_vec   .=  CS.Interp_Ct(WindFarm.u_ambient);  #Ct of each turbine
        CS.P_vec    .=  CS.Interp_P(WindFarm.u_ambient);   #P of each turbine
    elseif WindFarm.Turbine_Type=="IEA_15MW"
        WindFarm.D = 240;
        WindFarm.H = 150;
        data=load("04_Turbine_Data\\IEA_15MW.jld2")         # Load data from Input 
        WindFarm.Ct_Input;                                  # Thrustcoefficient vs. Windspeed
        WindFarm.P_Input;                                   # Power vs. Windspeed
        CS.Interp_Ct   =   LinearInterpolation(WindFarm.Ct_Input[:,1], WindFarm.Ct_Input[:,2]) #Interpolation function for Ct
        CS.Interp_P    =   LinearInterpolation(WindFarm.P_Input[:,1], WindFarm.P_Input[:,2]) #Interpolation function for P
        CS.Ct_vec   .=  CS.Interp_Ct(WindFarm.u_ambient);  #Ct of each turbine
        CS.P_vec    .=  CS.Interp_P(WindFarm.u_ambient);   #P of each turbine
    else
        error("ERROR: Wrong choice of turbine model in", WindFarm.Name,". Make sure to choose correctly between supported turbines.")
    end

    # If onluy two dimensional computation is conducted -> assign Z-Level to Hubheight
    # This bit needs some revision with respect to including different turbine (diameters) into the computation.
    if WindFarm.Dimensions == "2D"
        CS.XCoordinates .= CS.XCoordinates .* WindFarm.D
        CS.YCoordinates .= CS.YCoordinates .* WindFarm.D
        CS.ZCoordinates[1,:,1] .= WindFarm.H;
        CS.r        .= CS.YCoordinates  # Compute vector in radial & height direction for computation
    elseif WindFarm.Dimensions == "3D"
        CS.XCoordinates .= CS.XCoordinates .* WindFarm.D
        CS.YCoordinates .= CS.YCoordinates .* WindFarm.D
        CS.ZCoordinates[1,:,1] .= CS.ZCoordinates[1,:,1] .* WindFarm.D .+ WindFarm.H;
        CS.r            .= sqrt.(CS.YCoordinates.^2 .+ (CS.ZCoordinates.-WindFarm.H).^2) # Compute vector in radial & height direction for computation
    else
        error("EROOR: Wrong choice of dimensions in", WindFarm.Name,". Check Dimensions Input in the Input file.")
    end

end #LoadTurbineDATA

function LoadAtmosphericData!(WindFarm,CS)
#= This function loads all atmospheric data necessary for the computation
 It returns an updated WindFarm & CS struct with:
  1) Simple Computation: Wind & TI shear profile according to the height coordinates/ rotor resolution chosen by the user.
  2) AEP Computation: TBD
=#
 WindFarm.u_ambient_zprofile = zeros(1,CS.Real_Rotor_Res,1)#Assign right size to vector
 # Compute the ambient velocity log profile only for positive Z_Levels
 WindFarm.u_ambient_zprofile .= (CS.ZCoordinates .> 0) .* (WindFarm.u_ambient .* log.(CS.ZCoordinates ./ WindFarm.z_Surf) ./ log.(WindFarm.z_r / WindFarm.z_Surf))
 
 # Compute TI profile
    #TBDone!
end #LoadAtmosphericData

mutable struct ComputationStruct
    #Definition of struct (preassignment of arrays & space)
    #Coordinates/ Arrays
    XCoordinates::Array{Float64,3};
    YCoordinates::Array{Float64,3};
    ZCoordinates::Array{Float64,3};
    r::Array{Float64,3}; #Needed for single wake computation. Vector in radial & height direction.
    Real_Rotor_Res::Int; #Amount of computation points per rotor
    # Ambient data
    alpha_Comp::Float64;
    # Turbine specifics
    Yaw_Comp::Vector{Float64};      #Yawangle of each turbine
    Ct_vec::Array{Float64,3};       #Ct of each turbine
    P_vec::Array{Float64,3};        #P of each turbine 
    Interp_Ct::Any; #Interpolation function of thrustcoefficient
    Interp_P::Any;  #Interpolation function of power
    u_0_vec::Array{Float64,3};      #Inflow velocity of each turbine (Hubheight)
    u_0_vec_old::Array{Float64,3};  #Old Inflow velocity from last iteration 
    TI_0_vec::Array{Float64,3};  #Inflow turbulence intensity of each turbine (Hubheight)
    # Empirical values needed for Ishihara wake model
    k::Array{Float64,3};
    epsilon::Array{Float64,3};
    a::Array{Float64,3};
    b::Array{Float64,3};
    c::Array{Float64,3};
    d::Array{Float64,3};
    e::Array{Float64,3};
    f::Array{Float64,3};
    #Arrays needed for Wake computation
    sigma::Array{Float64,3};    #Wakewidth
    Delta_U::Array{Float64,3};  #Velocity deficit
    k1::Array{Float64,3};       #Parameter for turbulence computation
    k2::Array{Float64,3};       #Parameter for turbulence computation
    delta::Array{Float64,3};    #Parameter for turbulence computation
    Delta_TI::Array{Float64,3}; #Rotor-added turbulence
    Computation_Region_ID::Array{Bool,3}; #ID for limiting computation of the wake region
    #Arrays needed to superimpose
    U_Farm::Array{Float64,3};
    TI_Farm::Array{Float64,3};
    #Computation Parameters
    zeta::Float64 # termination criterion
end #mutable struct "ComputationStruct"

end #Module