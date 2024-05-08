# Module to initialise matrices and compute all necessary coordinates

module Initialisation_Module
using JLD2, Interpolations, LinearAlgebra#, MAT, LatinHypercubeSampling, PlotlyJS, Colors

export initCompArrays, LoadTurbineDATA, LoadAtmosphericData

function initCompArrays(WindFarm)
    println("Initialising arrays..")
    
    # Adjust user input for computation
    alpha_Comp  =   deg2rad(270 - WindFarm.alpha) #User Input logic: Geographical. Computation logic: Flow from left to right (270°(Western wind)==0°)
    Yaw_Comp    =   deg2rad.(270 .- WindFarm.Yaw)
    
    ### Compute coordinates (in D) for wake computation ###
    # Dimenisons: Relative X Coordinate, Relative Y Coordinate, Z Coordinate, Absolute turbine number
    XCoordinate = zeros(Float64, WindFarm.N , WindFarm.Y_Res, 1, WindFarm.N);  #Array for X coordinates of all points
    YCoordinate = zeros(Float64, WindFarm.N , WindFarm.Y_Res, 1, WindFarm.N);  #Array for Y Coordinates of all points
    Z_Levels    = zeros(Float64, 1, 1, WindFarm.Z_Res, 1);                    #Vector containing all height coordinates

    TMPYCoord   = zeros(Float64, WindFarm.Y_Res, WindFarm.N); #For yaw transformation of coordinates (Y coordinate)
    TMPXCoord   = zeros(Float64, WindFarm.Y_Res, WindFarm.N); #For yaw transformation of coordinates (X coordinate)

    #TMPZCoord   = zeros(Float64, WindFarm.RotorRes, 1, WindFarm.N); #For LHS distribution of Z coordinates

    ## Initialise temporary X & Y coordinate for yaw angle transformation
    #TMPYaw_X     = zeros(Float64, WindFarm.RotorRes, WindFarm.N);
    #TMPY_vector     = zeros(Float64, WindFarm.Y_Res);


    # Distribute the Y Coordinate 
    TMPY_vector     = LinRange(-1.5, 1.5, WindFarm.Y_Res) 

    
    for i in 1:WindFarm.N
    # Transform Y Coordinate according to yaw angle
    TMPYCoord[:,i]  = 0 .* sin(Yaw_Comp[i]) .+ TMPY_vector .* cos(Yaw_Comp[i]);
    TMPXCoord[:,i]  = 0 .* cos(Yaw_Comp[i]) .- TMPY_vector .* sin(Yaw_Comp[i]);
    end

    for i in 1:WindFarm.N
    # Create coordinate array of the structure: 1.Dim: Relative turbine, 2.Dim: RotorPoints, 3.Dim: Absolute turbine
    XCoordinate[:, :, 1, i] .= TMPXCoord' .+ WindFarm.x_vec .- WindFarm.x_vec[i];
    YCoordinate[:, :, 1, i]  = TMPYCoord' .+ WindFarm.y_vec .- WindFarm.y_vec[i];
    end
    Z_Levels[1,1,:,1] = LinRange(0, WindFarm.Z_Max, WindFarm.Z_Res);

    # Transform with respect to wind direction
    TMPX = zeros(Float64, WindFarm.N, WindFarm.Y_Res)
    for i in 1:WindFarm.N
        TMPX .= XCoordinate[:, :, 1, i];
        XCoordinate[:, :, 1, i] .=  TMPX .* cos(alpha_Comp) + YCoordinate[:, :, 1, i] .* sin(alpha_Comp);
        YCoordinate[:, :, 1, i] .= -TMPX .* sin(alpha_Comp) + YCoordinate[:, :, 1, i] .* cos(alpha_Comp);
    end
    
    # Transform with respect to yaw angle


#=COMMENTED OUT
3D Plot of Point!
    # Define colors
    colors = distinguishable_colors(WindFarm.N, colorant"blue")
    
    # Create an empty plot
    layout = Layout(
    scene=attr(
        xaxis=attr(
            nticks=4,
            range=[-20,20]
        ),
        yaxis=attr(
            nticks=4,
            range=[-20,20]
        ),
        zaxis=attr(
            nticks=4,
            range=[-1,7]
        ),
    ),
    scene_aspectratio=attr(x=2, y=2, z=0.2),
    width=700,
    margin=attr(
        r=20,
        l=10,
        b=10,
        t=10
    ),
    )
    
    scatter3d_plot = plot(layout)
    k=1;
    # Iterate over each row and add scatter traces to the plot
    for i in 1:WindFarm.N
        x = XCoordinate[i, :, 1, k]
        y = YCoordinate[i, :, 1, k]
        #z = ZCoordinate[i, :, 1, k]
        
        scatter_trace = scatter3d(x=x, y=y, z=z, mode="markers", marker_size=5, marker_color=colors[i], name="Set $i")
        add_trace!(scatter3d_plot, scatter_trace)
    end
    
    # Show plot
    display(scatter3d_plot)
#End Plot   =# 

    CS=ComputationStruct(   XCoordinate, YCoordinate, Z_Levels, zeros(WindFarm.N , WindFarm.Y_Res, WindFarm.Z_Res, WindFarm.N), alpha_Comp, Yaw_Comp,
                            zeros(1,1,1,WindFarm.N), zeros(1,1,1,WindFarm.N), (zeros(1,1,1,WindFarm.N) .+ WindFarm.u_ambient), 
                            zeros(1,1,1,WindFarm.N), zeros(1,1,1,WindFarm.N), zeros(1,1,1,WindFarm.N), zeros(1,1,1,WindFarm.N), zeros(1,1,1,WindFarm.N), 
                            zeros(1,1,1,WindFarm.N), zeros(1,1,1,WindFarm.N), zeros(1,1,1,WindFarm.N), zeros(size(XCoordinate)), zeros(size(XCoordinate)), 
                            zeros(size(XCoordinate)), zeros(size(XCoordinate)), zeros(size(XCoordinate)), zeros(size(XCoordinate))
                        )
 
    return WindFarm, CS

end #initCompArrays

function LoadTurbineDATA(WindFarm, CS)
#= This function loads all turbine data necessary for the computation
 It returns an update WindFarm & CS struct with Thrust coefficient and power curve, turbine diameter & hubheight. 
The ZCoordinate is also coorrected to have its origin at the Hubheigt of the turbine chosen for modelling. =#
    println("loading turbine data..")

    if WindFarm.VestasV80==true && WindFarm.NREL_5MW==false
        WindFarm.D = 80;
        WindFarm.H = 70;
        data = load("04_Turbine_Data\\VestasV80_2MW.jld2")  # Load data from Input 
        WindFarm.Ct_Input = data["CT_Input"];               # Thrustcoefficient vs. Windspeed
        WindFarm.P_Input = data["P_Input"];                 # Power vs. Windspeed
        Interp_Ct   =   LinearInterpolation(WindFarm.Ct_Input[:,1], WindFarm.Ct_Input[:,2]) #Interpolation function for Ct
        Interp_P    =   LinearInterpolation(WindFarm.P_Input[:,1], WindFarm.P_Input[:,2]) #Interpolation function for P
        CS.Ct_vec   .=  Interp_Ct(WindFarm.u_ambient);  #Ct of each turbine
        CS.P_vec    .=  Interp_P(WindFarm.u_ambient);   #P of each turbine 

        CS.Z_Levels .= CS.Z_Levels ./WindFarm.D;  #Adjust ZCoodinates for computation according to Hubheight
    elseif WindFarm.NREL_5MW==true && WindFarm.VestasV80==false
        WindFarm.D = 126;
        WindFarm.H = 90;
        data=load("04_Turbine_Data\\NREL_5MW.jld2")         # Load data from Input 
        WindFarm.Ct_Input;                                  # Thrustcoefficient vs. Windspeed
        WindFarm.P_Input;                                   # Power vs. Windspeed
        Interp_Ct   =   LinearInterpolation(WindFarm.Ct_Input[:,1], WindFarm.Ct_Input[:,2]) #Interpolation function for Ct
        Interp_P    =   LinearInterpolation(WindFarm.P_Input[:,1], WindFarm.P_Input[:,2]) #Interpolation function for P
        CS.Ct_vec   .=  Interp_Ct(WindFarm.u_ambient);  #Ct of each turbine
        CS.P_vec    .=  Interp_P(WindFarm.u_ambient);   #P of each turbine
        
        CS.Z_Levels .= CS.Z_Levels .+ WindFarm.H./WindFarm.D;  #Adjust ZCoodinates for computation according to Hubheight 
    else
        error("ERROR: Wrong choice of turbine model in", WindFarm.Name,"Make sure to choose one but not more.")
    end

    return WindFarm, CS
end #LoadTurbineDATA

function LoadAtmosphericData(WindFarm,CS)
#= This function loads all atmospheric data necessary for the computation
 It returns an updated WindFarm & CS struct with:
  1) Simple Computation: Wind & TI shear profile according to the height coordinates/ rotor resolution chosen by the user.
  2) AEP Computation: TBD
=#
 println("loading atmospheric data..")
 WindFarm.u_ambient_zprofile=Array{Float64,4}(undef,1,1,WindFarm.Z_Res,1) #Assign right size to vector

 # Compute ambient velocity for each coordinate point according to wind shear profile
 WindFarm.u_ambient_zprofile .= ifelse.(CS.Z_Levels.>0, WindFarm.u_ambient .* log.(CS.Z_Levels.*WindFarm.D./WindFarm.z_Surf)./log.(WindFarm.z_r./WindFarm.z_Surf),0);

 # Compute TI profile
    #TBDone!

 return WindFarm, CS
end #LoadAtmosphericData

mutable struct ComputationStruct
    #Definition of struct (preassignment of arrays & space)
    #Coordinates/ Arrays
    XCoordinates::Array{Float64,4};
    YCoordinates::Array{Float64,4};
    Z_Levels::Array{Float64,4};
    r::Array{Float64,4}; #Needed for single wake computation. Vector in radial & height direction.
    # Ambient data
    alpha_Comp::Float64;
    # Turbine specifics
    Yaw_Comp::Vector{Float64};  #Yawangle of each turbine
    Ct_vec::Array{Float64,4};    #Ct of each turbine
    P_vec::Array{Float64,4};     #P of each turbine 
    u_0_vec::Array{Float64,4};   #Inflow velocity of each turbine  
    # Empirical values needed for Ishihara wake model
    k::Array{Float64,4};
    epsilon::Array{Float64,4};
    a::Array{Float64,4};
    b::Array{Float64,4};
    c::Array{Float64,4};
    d::Array{Float64,4};
    e::Array{Float64,4};
    f::Array{Float64,4};
    #Arrays needed for Wake computation
    sigma::Array{Float64,4};    #Wakewidth
    Delta_U::Array{Float64,4};  #Velocity deficit
    k1::Array{Float64,4};       #Parameter for turbulence computation
    k2::Array{Float64,4};       #Parameter for turbulence computation
    delta::Array{Float64,4};    #Parameter for turbulence computation
    Delta_TI::Array{Float64,4}; #Rotor-added turbulence

end #mutable struct "ComputationStruct"

end #Module