# Module to initialise matrices and compute all necessary coordinates

module Initialisation_Module
using LatinHypercubeSampling#, PlotlyJS, Colors

export initCompArrays

function initCompArrays(WindFarm)
    println("Initialising arrays..")
    
    # Adjust user input for computation
    alpha_Comp  =   deg2rad(270 - WindFarm.alpha) #User Input logic: Geographical. Computation logic: Flow from left to right (270°(Western wind)==0°)
    Yaw_Comp    =   deg2rad.(270 .- WindFarm.Yaw)

    ### Compute coordinates (in D) for wake computation ###
    # Dimenisons: Turbine's relative coordinate, Number of Points, Turbines number
    XCoordinate = zeros(Float64, WindFarm.N-1 , WindFarm.RotorRes, WindFarm.N);
    YCoordinate = zeros(Float64, WindFarm.N-1 , WindFarm.RotorRes, WindFarm.N);
    ZCoordinate = zeros(Float64, WindFarm.N-1 , WindFarm.RotorRes, WindFarm.N) .+ (WindFarm.H/WindFarm.D); #Preassign hubheight to all coordinates (Normalised to D)
    TMPYCoord   = zeros(Float64, WindFarm.RotorRes, WindFarm.N); #For LHS distribution of Y coordinates
    TMPZCoord   = zeros(Float64, WindFarm.RotorRes, WindFarm.N); #For LHS distribution of Z coordinates
    TMPRCoord   = zeros(Float64, WindFarm.RotorRes);    #For LHS distribution of Y&Z (First computation in polars then transformation to cartesian)
    TMPphiCoord = zeros(Float64, WindFarm.RotorRes);    #For LHS distribution of Y&Z (First computation in polars then transformation to cartesian)
    LHCPlan     = zeros(Float64, WindFarm.RotorRes,2);  #For LHS distribution of Y&Z

    ## Initialise temporary X & Y coordinate for yaw angle transformation
    TMPYaw_X     = zeros(Float64, WindFarm.RotorRes, WindFarm.N);
    TMPYaw_Y     = zeros(Float64, WindFarm.RotorRes);

    ## Determine LHS points of each turbines Y & Z Coordinate (circular)
    #LHS for spanwise & height distribution of points (y&z coordinates) are first computed as polar coordinates.
    for i in 1:WindFarm.N
    LHCPlan = randomLHC(WindFarm.RotorRes,2)./WindFarm.RotorRes; #LHC plan for polars (r & phi).
    TMPRCoord   = LHCPlan[:,1].*0.5;      #Radius always [0:0.5] (coordinates normalised to D)
    TMPphiCoord = LHCPlan[:,2].*2pi;      #Angle in RAD
    
    # Now convert & store coordinates as cartesian
    TMPYCoord[:,i] = TMPRCoord .* cos.(TMPphiCoord);
    TMPZCoord[:,i] = TMPRCoord .* sin.(TMPphiCoord) .+ WindFarm.H/WindFarm.D;

    # Transform Y Coordinate according to yaw angle
    TMPYaw_Y = TMPYCoord[:,i];
    TMPYCoord[:,i] = 0 .* sin(Yaw_Comp[i]) .+ TMPYCoord[:,i] .* cos(Yaw_Comp[i]);
    TMPYaw_X[:,i]  = 0 .* cos(Yaw_Comp[i]) .- TMPYaw_Y       .* sin(Yaw_Comp[i]);
    end

    ## Find coordinates of all turbines relative to each other
    XCoordinate[:, :, 1] .= TMPYaw_X[:,2:end]' .+ WindFarm.x_vec[2:end] .- WindFarm.x_vec[1];
    YCoordinate[:, :, 1]  = TMPYCoord[:,2:end]' .+ WindFarm.y_vec[2:end] .- WindFarm.y_vec[1];
    ZCoordinate[:, :, 1]  = TMPZCoord[:,2:end]';

    # Create coordinate array of the structure: 1.Dim: Relative turbine, 2.Dim: RotorPoints, 3.Dim: Absolute turbine
    for i in 2:WindFarm.N
        XCoordinate[:, :, i] .= vcat(TMPYaw_X[:,1:i-1]', TMPYaw_X[:,i+1:end]') .+ vcat(WindFarm.x_vec[1:i-1], WindFarm.x_vec[i+1:end]) .- WindFarm.x_vec[i];
        YCoordinate[:, :, i] .= vcat(TMPYCoord[:,1:i-1]', TMPYCoord[:,i+1:end]') .+ vcat(WindFarm.y_vec[1:i-1], WindFarm.y_vec[i+1:end]) .- WindFarm.y_vec[i];
        ZCoordinate[:, :, i] .= vcat(TMPZCoord[:,1:i-1]', TMPZCoord[:,i+1:end]');
    end
    
    # Transform with respect to wind direction
    TMPX = zeros(Float64, WindFarm.N-1, WindFarm.RotorRes)
    for i in 1:WindFarm.N
        TMPX .= XCoordinate[:, :, i];
        XCoordinate[:, :, i] .=  TMPX .* cos(alpha_Comp) + YCoordinate[:, :, i] .* sin(alpha_Comp);
        YCoordinate[:, :, i] .= -TMPX .* sin(alpha_Comp) + YCoordinate[:, :, i] .* cos(alpha_Comp);
    end
    
    # Transform with respect to yaw angle


#= COMMENTED OUT
3D Plot of Point!
    # Define colors
    colors = distinguishable_colors(WindFarm.N-1, colorant"blue")
    
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
    for i in 1:WindFarm.N-1
        x = XCoordinate[i, :,k]
        y = YCoordinate[i, :,k]
        z = ZCoordinate[i, :,k]
        
        scatter_trace = scatter3d(x=x, y=y, z=z, mode="markers", marker_size=5, marker_color=colors[i], name="Set $i")
        add_trace!(scatter3d_plot, scatter_trace)
    end
    
    # Show plot
    display(scatter3d_plot)
#End Plot   =# 


    CD=ComputationStruct(XCoordinate, YCoordinate, ZCoordinate, alpha_Comp, Yaw_Comp)
    return CD
end #initCompArrays


mutable struct ComputationStruct
    #Definition of struct (preassignment of arrays & space)
    XCoodinates::Array{Float64,3};
    YCoordinates::Array{Float64,3};
    ZCoordinates::Array{Float64,3};
    alpha_Comp::Float64;
    Yaw_Comp::Vector{Float64};
end #mutable struct "ComputationStruct"

end #Module