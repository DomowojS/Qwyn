# Module to initialise matrices and compute all necessary coordinates

module Initialisation_Module
using JLD2, Interpolations, LinearAlgebra, Random, StatsBase, Statistics, FastGaussQuadrature

export initCompArrays, LoadTurbineDATA!, LoadAtmosphericData, FindStreamwiseOrder!, initGraphicArrays, correctGS

function initCompArrays(WindFarm)
# Initialises/ preallocates all arrays needed for computation
    ### XYZ initialisation    
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
            Y_vec, Z_vec, W_vec = generate_rotor_grid(WindFarm.Rotor_Res)
        elseif WindFarm.Rotor_Discretization == "fibonacci"
            #Generates a fibonacci-lattice distributed grid. Accurate representation for significantly lower amount of points.
            Y_vec, Z_vec, W_vec = generate_fibonacci(WindFarm.Rotor_Res)
        elseif WindFarm.Rotor_Discretization == "smart_grid"
            #Generates circles with even distribution
            Y_vec, Z_vec, W_vec = generate_smart_rotorgrid(WindFarm.Rotor_Res)
        else
            error("ERROR: Wrong choice of rotor discretization function for", WindFarm.Name,". Check variable Rotor_Discretization. Currently allowed entries: gridded, fibonacci, smart_grid.")
        end

        # Find and stor actual amount of points per rotor (after grid generation)
        Real_Rotor_Res = length(Y_vec);

        ### Compute coordinates (in D) for wake computation ###
        # Dimenisons: Relative X Coordinate, Relative Y Coordinate, Z Coordinate, Absolute turbine number
        XCoordinate = zeros(Float64, WindFarm.N , 1, WindFarm.N);  #Array for X coordinates of all points
        YCoordinate = zeros(Float64, WindFarm.N , Real_Rotor_Res, WindFarm.N);  #Array for Y Coordinates of all points
        ZCoordinate = zeros(Float64, 1 , Real_Rotor_Res, 1);  #Vector containing all height coordinates
        RotorPointWeights= zeros(Float64, 1 , Real_Rotor_Res, 1);  #Vector containing all Rotor point weights

        for i in 1:WindFarm.N
        # Create coordinate array of the structure: 1.Dim: Relative turbine, 2.Dim: RotorPoints, 3.Dim: Absolute turbine
            XCoordinate[:, :, i] .= XCoordinate[:, :, i] .+ WindFarm.x_vec .- WindFarm.x_vec[i];
            for j in 1:WindFarm.N
                YCoordinate[j, :, i] .= Y_vec .+ WindFarm.y_vec[j];
            end
            YCoordinate[:,:,i] .= YCoordinate[:,:,i] .- WindFarm.y_vec[i]
        end
        ZCoordinate[1,:,1] .= Z_vec
        RotorPointWeights[1,:,1] .= W_vec

    ### Generate grid for convection velocity computation (needed for momentum conserving superposition)
        if WindFarm.Superpos == "Momentum_Conserving"
        # generate basic grid vector 

        if WindFarm.Uc_Res < 4; WindFarm.Uc_Res=4, println("Uc_Res corrected to minimal value 4") end #correction if resolution is chosen too small

        Y_Uc, Z_Uc, WindFarm.Uc_Res = generate_grid_for_Uc(minimum(WindFarm.y_vec)-1, maximum(WindFarm.y_vec)+1, 0.1, 2.0, WindFarm.Uc_Res) #Height normalised -> to be scaled after turbine diameter and hub height read in in "LoadTurbineData"
        
        #Preassign arrays for coordinates (the x coordinates can be reused from the rotor points array)
        Y_for_Uc            =   zeros(Float64, WindFarm.N , WindFarm.Uc_Res, WindFarm.N);
        Z_for_Uc            =   zeros(Float64, 1 , WindFarm.Uc_Res, 1);
        r_for_Uc            =   zeros(WindFarm.N, WindFarm.Uc_Res, WindFarm.N);
        #Preassign arrays for result
        
        
        #Create coordinate arrays
        for ii in 1:WindFarm.N
            # Create coordinate array of the structure: 1.Dim: Relative turbine, 2.Dim: RotorPoints, 3.Dim: Absolute turbine
                for jj in 1:WindFarm.N
                    Y_for_Uc[jj, :, ii] .= Y_Uc .+ WindFarm.y_vec[jj];
                end
                Y_for_Uc[:,:,ii] .= Y_for_Uc[:,:,ii] .- WindFarm.y_vec[ii]
            end
            Z_for_Uc[1,:,1] .= Z_Uc
        
        #Preassignment of additional arrays needed for computation
        u_ambient_for_Uc            = zeros(1, WindFarm.Uc_Res, 1)
        Delta_U_for_Uc              = zeros(WindFarm.N, WindFarm.Uc_Res, WindFarm.N)
        Mixed_wake_for_Uc           = zeros(WindFarm.N, WindFarm.Uc_Res, 1)
        Comutation_Region_ID_for_Uc = trues(WindFarm.N, WindFarm.Uc_Res,WindFarm.N)
        sigma_for_Uc                = zeros(WindFarm.N, WindFarm.Uc_Res, WindFarm.N)
        sigma_m_for_Uc              = zeros(WindFarm.N, WindFarm.Uc_Res, WindFarm.N)
        Lambda_for_Uc               = zeros(1, WindFarm.Uc_Res, WindFarm.N)
        k1_for_Uc                   = zeros(WindFarm.N, WindFarm.Uc_Res, WindFarm.N)  
        k2_for_Uc                   = zeros(WindFarm.N, WindFarm.Uc_Res, WindFarm.N)  
        delta_for_Uc                = zeros(1, WindFarm.Uc_Res, 1)  
        Delta_TI_for_Uc             = zeros(WindFarm.N, WindFarm.Uc_Res, WindFarm.N)
        weighting_Factor_for_Uc     = zeros(WindFarm.N, WindFarm.Uc_Res, WindFarm.N)
        
        else

        #Dummies for struct definition in case other superposition is used 
        Y_for_Uc                    =   zeros(1,1,1)
        Z_for_Uc                    =   zeros(1,1,1) 
        r_for_Uc                    =   zeros(1,1,1)
        u_ambient_for_Uc            =   zeros(1,1,1)
        Delta_U_for_Uc              =   zeros(1,1,1)
        Mixed_wake_for_Uc           =   zeros(1,1,1)
        Comutation_Region_ID_for_Uc =   trues(1,1,1)
        sigma_for_Uc                =   zeros(1,1,1)
        sigma_m_for_Uc              =   zeros(1,1,1)
        Lambda_for_Uc               =   zeros(1,1,1)
        k1_for_Uc                   =   zeros(1,1,1)
        k2_for_Uc                   =   zeros(1,1,1)
        delta_for_Uc                =   zeros(1,1,1)
        Delta_TI_for_Uc             =   zeros(1,1,1)
        weighting_Factor_for_Uc     =   zeros(1,1,1)
        end    
        
    #= Transform with respect to yaw angle
    TBDone!!
    =#


    # Create struct which holds all computation arrays
    CS=ComputationStruct(WindFarm.name, 0.0, XCoordinate, YCoordinate, ZCoordinate, RotorPointWeights, zeros(WindFarm.N , Real_Rotor_Res, WindFarm.N), Real_Rotor_Res,  alpha_Comp, Yaw_Comp,
                            zeros(1,1,WindFarm.N), zeros(1,1,WindFarm.N), 0, 0, 0, (zeros(1,1,WindFarm.N) .+ WindFarm.u_ambient), (zeros(1,1,WindFarm.N) .+ WindFarm.TI_a),
                            zeros(1,1,WindFarm.N), zeros(1,1,WindFarm.N), zeros(1,1,WindFarm.N), zeros(1,1,WindFarm.N), zeros(1,1,WindFarm.N), zeros(1,1,WindFarm.N), zeros(1,1,WindFarm.N), zeros(1,1,WindFarm.N), 
                            zeros(WindFarm.N,Real_Rotor_Res,WindFarm.N), zeros(WindFarm.N,Real_Rotor_Res,WindFarm.N), zeros(WindFarm.N,Real_Rotor_Res,WindFarm.N), zeros(WindFarm.N,Real_Rotor_Res,WindFarm.N), 
                            zeros(1,Real_Rotor_Res,1), zeros(WindFarm.N,Real_Rotor_Res,WindFarm.N), trues(WindFarm.N,Real_Rotor_Res,WindFarm.N), zeros(WindFarm.N,1,WindFarm.N), zeros(WindFarm.N,1,1), 
                            zeros(WindFarm.N,1,1), zeros(WindFarm.N,Real_Rotor_Res,WindFarm.N), Delta_U_for_Uc, Mixed_wake_for_Uc, Y_for_Uc, Z_for_Uc, r_for_Uc, u_ambient_for_Uc, Comutation_Region_ID_for_Uc, 
                            sigma_for_Uc, sigma_m_for_Uc, Lambda_for_Uc, k1_for_Uc, k2_for_Uc, delta_for_Uc, Delta_TI_for_Uc, weighting_Factor_for_Uc, zeros(1,1,WindFarm.N), zeros(1,Real_Rotor_Res,WindFarm.N), 
                            zeros(WindFarm.N,Real_Rotor_Res,WindFarm.N), zeros(WindFarm.N,Real_Rotor_Res,1), zeros(WindFarm.N,Real_Rotor_Res,1), trues(WindFarm.N), 
                            trues(WindFarm.N), falses(WindFarm.N),0, zeros(WindFarm.N), zeros(WindFarm.N), zeros(WindFarm.N,1,WindFarm.N)
                        )
 
    return WindFarm, CS

end #initCompArrays

function LoadTurbineDATA!(WindFarm, CS)
#= This function loads all turbine data necessary for the computation
 It returns an update WindFarm & CS struct with Thrust coefficient and power curve, turbine diameter & hubheight. 
The ZCoordinate is also coorrected to have its origin at the Hubheigt of the turbine chosen for modelling. =#
    if WindFarm.Turbine_Type=="VestasV80"
        WindFarm.D = 80;
        WindFarm.H = 70;
        data = load("04_Turbine_Data/VestasV80_2MW.jld2")  # Load data from Input 
        WindFarm.Ct_Input = data["CT_Input"];               # Thrustcoefficient vs. Windspeed
        WindFarm.P_Input  = data["P_Input"];                 # Power vs. Windspeed
        
        # Correct input values in lower and upper region
        lower_correctionCT  = [0 0; 3.99 0];
        upper_correctionCT  = [25.0135 0; 50 0];
        WindFarm.Ct_Input = vcat(lower_correctionCT, WindFarm.Ct_Input)
        WindFarm.Ct_Input = vcat(WindFarm.Ct_Input, upper_correctionCT)

        lower_correctionP  = [0 0; 3.98 0];
        upper_correctionP  = [25.0135 0; 50 0];
        WindFarm.P_Input = vcat(lower_correctionP, WindFarm.P_Input)
        WindFarm.P_Input = vcat(WindFarm.P_Input, upper_correctionP)

    elseif WindFarm.Turbine_Type=="NREL_5MW"
        WindFarm.D = 126;
        WindFarm.H = 90;
        data=load("04_Turbine_Data/NREL_5MW.jld2")         # Load data from Input 
        WindFarm.Ct_Input;                                  # Thrustcoefficient vs. Windspeed
        WindFarm.P_Input;                                   # Power vs. Windspeed
    elseif WindFarm.Turbine_Type=="DTU_10MW"
        WindFarm.D = 119;
        WindFarm.H = 178.3;
        data=load("04_Turbine_Data/DTU_10MW.jld2")         # Load data from Input 
        WindFarm.Ct_Input;                                  # Thrustcoefficient vs. Windspeed
        WindFarm.P_Input;                                   # Power vs. Windspeed
    elseif WindFarm.Turbine_Type=="IEA_15MW"
        WindFarm.D = 240;
        WindFarm.H = 150;
        data=load("04_Turbine_Data/IEA_15MW.jld2")         # Load data from Input 
        WindFarm.Ct_Input;                                  # Thrustcoefficient vs. Windspeed
        WindFarm.P_Input;                                   # Power vs. Windspeed
    else
        error("ERROR: Wrong choice of turbine model in", WindFarm.Name,". Make sure to choose correctly between supported turbines.")
    end
    #Define Interpolation functions and preassign initial Ct and P values
    CS.Interp_Ct   =   LinearInterpolation(WindFarm.Ct_Input[:,1], WindFarm.Ct_Input[:,2]) #Interpolation function for Ct
    CS.Interp_P    =   LinearInterpolation(WindFarm.P_Input[:,1], WindFarm.P_Input[:,2]) #Interpolation function for P
    CS.Ct_vec   .=  CS.Interp_Ct(WindFarm.u_ambient);  #Ct of each turbine
    CS.P_vec    .=  CS.Interp_P(WindFarm.u_ambient);   #P of each turbine 

        ### Correction of XYZ to rotor diameter D and Hub height
        # This bit needs some revision with respect to including different turbine (diameters) into the computation.
            # Rotor points
            CS.XCoordinates         .= CS.XCoordinates .* WindFarm.D
            CS.YCoordinates         .= CS.YCoordinates .* WindFarm.D
            CS.ZCoordinates         .= CS.ZCoordinates .* WindFarm.D .+ WindFarm.H;
            CS.r                    .= sqrt.(CS.YCoordinates.^2 .+ (CS.ZCoordinates.-WindFarm.H).^2) # Compute vector in radial & height direction for computation
            # Points for convection velocity
            if WindFarm.Superpos == "Momentum_Conserving"
                CS.Y_for_Uc         .= CS.Y_for_Uc .* WindFarm.D
                CS.Z_for_Uc         .= CS.Z_for_Uc .* WindFarm.D;
                CS.r_for_Uc         .= sqrt.(CS.Y_for_Uc.^2 .+ (CS.Z_for_Uc .- WindFarm.H).^2) # Compute vector in radial & height direction for computation
            end

end #LoadTurbineDATA

function LoadAtmosphericData(WindFarm, CS)
#= This function loads all atmospheric data necessary for the computation
 It returns an updated WindFarm & CS struct with:
  1) Simple Computation: Wind & TI shear profile according to the height coordinates/ rotor resolution chosen by the user.
  2) AEP Computation: TBD
=#
    u_ambient_zprofile  = zeros(1, CS.Real_Rotor_Res,1)#Assign right size to vector
    # Compute the ambient velocity log profile only for positive Z_Levels
    u_ambient_zprofile .= (CS.ZCoordinates .> 0) .* (WindFarm.u_ambient .* log.(CS.ZCoordinates ./ WindFarm.z_Surf) ./ log.(WindFarm.z_r / WindFarm.z_Surf))

    # Initial values
    CS.u_0_vec .= mean(u_ambient_zprofile) #Inflow velocity of each turbine

    # Compute the ambient velocity also for points needed for Uc
    if WindFarm.Superpos == "Momentum_Conserving"
        CS.u_ambient_for_Uc         .= (CS.Z_for_Uc .> 0) .* (WindFarm.u_ambient .* log.(CS.Z_for_Uc ./ WindFarm.z_Surf) ./ log.(WindFarm.z_r / WindFarm.z_Surf))
    end

 # Correct inflow velocity


 # Compute TI profile
    #TBDone! -> According to Demetri (Bouris-NTUA) - not neccesarily needed.

    return u_ambient_zprofile, CS

end #LoadAtmosphericData

function FindStreamwiseOrder!(WindFarm, CS)
#= This function evaluates the wind turbine's coordinates in streamwise direction
and assigns a computation order =#

    CS.StreamwiseOrder    =   sortperm(WindFarm.x_vec); # Find streamwise order
    j=1; #Priority level counter
    ii=1;#Counter for turbines within level
    
    for i = 1:WindFarm.N

        if i > 1 # 2nd to nth turbine 
            
            if any(WindFarm.y_vec[CS.StreamwiseOrder[i]] .< (2 .* WindFarm.x_vec[CS.StreamwiseOrder[i]] .+ (WindFarm.y_vec[CS.StreamwiseOrder[ii:i-1]].- 2 .* WindFarm.x_vec[CS.StreamwiseOrder[ii:i-1]]))) == true && 
                any(WindFarm.y_vec[CS.StreamwiseOrder[i]] .> (-2 .* WindFarm.x_vec[CS.StreamwiseOrder[i]] .+ (WindFarm.y_vec[CS.StreamwiseOrder[ii:i-1]].+ 2 .* WindFarm.x_vec[CS.StreamwiseOrder[ii:i-1]]))) == true &&
                any(abs.(WindFarm.y_vec[CS.StreamwiseOrder[i]] .- WindFarm.y_vec[CS.StreamwiseOrder[ii:i-1]]) .< 2) == true #Check if new turbine is in shadowing region of any previous turbine
                ii=i;
                j = j+1; #Raise counter to one - new, lower priority level.

            end
        
        end

    CS.CompOrder[CS.StreamwiseOrder[i]]=j; #Assign order by streamwise position

    end#for
end#FindStreamwiseOrder

### Subfunctions & Structdefinitions #####
function generate_rotor_grid(totalPoints::Int)
# This function generates a grid of points in the Y-Z plane to represent the turbine's rotor's.
    # Handle special cases for 1 to 4 points
    if totalPoints == 1
        return [0.0], [0.0]
    elseif totalPoints == 2
        angle = π / 2
        return [0.5*cos(angle), -0.5*cos(angle)], [0.5*sin(angle), -0.5*sin(angle)]
    elseif totalPoints == 3
        angles = range(0, stop=2*π, length=4)[1:3]
        return 0.5*cos.(angles), 0.5*sin.(angles)
    elseif totalPoints == 4
        angles = range(π/4, stop=2*π, length=5)[1:4]
        return 0.5*cos.(angles), 0.5*sin.(angles)
    end

    # General case for totalPoints >= 5
    radius = 0.5
    numPointsPerAxis = ceil(Int, sqrt(4 * totalPoints / π))

    # Create a sufficiently dense grid
    x = range(-radius, stop=radius, length=numPointsPerAxis)
    y = range(-radius, stop=radius, length=numPointsPerAxis)
    X, Y = [xi for xi in x, yi in y], [yi for xi in x, yi in y]

    # Flatten the grid matrices
    X = vec(X)
    Y = vec(Y)

    # Calculate distance from the origin to each point
    distances = sqrt.(X.^2 .+ Y.^2)

    # Filter out points that are outside the rotor area of radius 0.5
    insiderotor = distances .<= radius

    # Keep only the points inside the circle
    Y = Y[insiderotor]
    Z = X[insiderotor]

    # Adjust if we have fewer points than needed
    while length(Y) < totalPoints
        # Increase the grid density
        numPointsPerAxis += 1
        x = range(-radius, stop=radius, length=numPointsPerAxis)
        y = range(-radius, stop=radius, length=numPointsPerAxis)
        X, Y = [xi for xi in x, yi in y], [yi for xi in x, yi in y]
        X = vec(X)
        Y = vec(Y)
        distances = sqrt.(X.^2 .+ Y.^2)
        insiderotor = distances .<= radius
        Y = Y[insiderotor]
        Z = X[insiderotor]
    end

    # If we have more points than needed, sample without replacement to get the exact number
    if length(Y) > totalPoints
        indices = sample(1:length(Y), totalPoints, replace=false)
        Y = Y[indices]
        Z = Z[indices]
    end

    # Even weight distribution
    W = fill(1 / totalPoints, totalPoints)  # Equal weighting for Fibonacci points

    return Y, Z, W
end # generate_rotor_grid
  
function generate_fibonacci(totalPoints::Int)
# This function generates a fibonacci-lattice grid of points int the Y-Z plane to represent the turbine's rotor's.
# Golden ratio
    φ = (sqrt(5) + 1) / 2

    # Generate points using the Fibonacci lattice method
    indices = 1:totalPoints

    # Calculate radii and angles
    r = sqrt.(indices / totalPoints) * 0.5  # radius scaled to the circle with radius 0.5
    θ = 2 * π * φ * indices

    # Convert polar coordinates to Cartesian coordinates
    Y = r .* cos.(θ)
    Z = r .* sin.(θ)

    # Even weight distribution
    W = fill(1 / totalPoints, totalPoints)  # Equal weighting for Fibonacci points
    
    return Y, Z, W
end #generate_fibonacci

function generate_smart_rotorgrid(n::Int)
#=
   Generate integration points for rotor averaging.
    
   Parameters:
   - n: Number of points. Allowed values: 1, 4, 6, 7, 9, 12, 21, or any value > 21.
   
   Returns:
   - Y, Z: Cartesian coordinates of the integration points.
   - weights: Corresponding integration weights.
   
   Throws an error if an unsupported number of points is provided.
=#
    
    if n ∈ [4, 7, 9, 21]
    # Predefined quadrature points and weights for CGI
    pm = [[-0.5, -0.5], [-0.5, 0.5], [0.5, -0.5], [0.5, 0.5]]  # Base square points
    YZ = Vector{Tuple{Float64, Float64}}()  # Ensure it can store floating-point values
    end

    if n == 1 #One center point at the hub
        return [0.0], [0.0], [1.0]  # Single point with full weight
    elseif n == 4 #4 points distributed using the circular gauss integration scheme
        append!(YZ, [(0.5 * p[1], 0.5 * p[2]) for p in pm])
        W = fill(1 / 4, 4)
        return first.(YZ), last.(YZ), W

    elseif n == 6 #6 points evenly distributed at r=2/3R (normalized to D)
        theta = range(-π, π; length=7)[1:6]  # 6 evenly spaced angles
        r = 2 / 3 .*0.5
        Y = r * cos.(-theta .- π / 2)
        Z = r * sin.(-theta .- π / 2)
        W = fill(1 / 6, 6)  # Equal weighting
        return Y, Z, W

    elseif n == 7 #7 points distributed using the circular gauss integration scheme
        push!(YZ, (0.0, 0.0))
        append!(YZ, [(-sqrt(2 / 3).*0.5, 0.0), (sqrt(2 / 3).*0.5, 0.0)]) # Normalized to D
        append!(YZ, [(sqrt(1 / 6) * p[1], sqrt(1 / 2) * p[2]) for p in pm])
        W = [1 / 4; fill(1 / 8, 6)]
        return first.(YZ), last.(YZ), W

    elseif n == 9 #9 points distributed using the circular gauss integration scheme
        push!(YZ, (0.0, 0.0))
        append!(YZ, [(-0.5, 0.0), (0.5, 0.0), (0.0, -0.5), (0.0, 0.5)])
        append!(YZ, [(0.5 * p[1], 0.5 * p[2]) for p in pm])
        W = [1 / 6, 1 / 24, 1 / 24, 1 / 24, 1 / 24, 1 / 6, 1 / 6, 1 / 6, 1 / 6]
        return first.(YZ), last.(YZ), W

    elseif n == 12 #12 point distributed using the gauss quadratute 

        # Gauss-Legendre quadrature nodes and weights for [-1,1] interval
        y, weights_y = gausslegendre(4)
        z, weights_z = gausslegendre(4)

        # Create a 2D meshgrid of points
        Y = repeat(y, inner=4)
        Z = repeat(z, outer=4)
        W = repeat(weights_y, inner=4) .* repeat(weights_z, outer=4)  # Compute 2D weights

        # Keep only points inside the rotor
        mask = Y.^2 + Z.^2 .< 1      
        Y, Z, W = Y[mask].*0.5, Z[mask].*0.5, W[mask] ./ sum(W[mask])  # Apply mask correctly

        return Y, Z, W
    
    elseif n == 21 #21 points distributed using the circular gauss integration scheme
        θ_1 = [2 * π * k / 10 for k in 1:10]
        θ_2 = [2 * π * k / 10 for k in 1:10]
        r_1 = sqrt.((6 - sqrt(6)) / 10) .*0.5
        r_2 = sqrt.((6 + sqrt(6)) / 10) .*0.5
        
        Y = vcat(0.0, r_1 .* cos.(θ_1), r_2 .* cos.(θ_2))
        Z = vcat(0.0, r_1 .* sin.(θ_1), r_2 .* sin.(θ_2))
        W = vcat(1 / 9, fill((16 + sqrt(6)) / 360, 10), fill((16 - sqrt(6)) / 360, 10))
        return Y, Z, W

    elseif n > 21 #>21 points distributed using the fibonacci lattice method
        Y, Z, W = generate_fibonacci(n)
        return Y, Z, W

    else
        error("Invalid number of Rotorpoints! For Smart grid generation choose Rotor_Res = 1, 4, 6, 7, 8, 9, 21, or n > 21.")
    end

end

function generate_grid_for_Uc(min_y::Float64, max_y::Float64, min_z::Float64, max_z::Float64, num_points::Int)
#Grid points for Uc
    if num_points < 4
        error("Number of points must be at least 4 to ensure points at each corner of the rectangle.")
    end

    # Calculate the number of points along each axis to form a grid
    num_points_side = ceil(Int, sqrt(num_points))

    # Generate evenly spaced points for Y and Z coordinates, including the maximum values
    y_coords = range(min_y, stop=max_y, length=num_points_side)
    z_coords = range(min_z, stop=max_z, length=num_points_side)

    # Create meshgrid
    Y = repeat(y_coords, inner=num_points_side)
    Z = repeat(z_coords, outer=num_points_side)

    return Y, Z, length(Y)
end

function initGraphicArrays(WindFarm)
# Initialises/ preallocates all arrays needed for Graphikcomputation
    # Adjust user input for computation
    alpha_Comp  =   deg2rad(270 - WindFarm.alpha) #User Input logic: Geographical. Computation logic: Flow from left to right (270°(Western wind)==0°)
    Yaw_Comp    =   deg2rad.(270 .- WindFarm.Yaw)

    # Create primal grids for x y and z coordinates
    X_vec, Y_vec, Z_vec, Real_Rotor_Res = generate_graphic_grid(WindFarm)
    

        ### Compute coordinates (in D) for wake computation ###
        # Dimenisons: Relative X Coordinate, Relative Y Coordinate, Z Coordinate, Absolute turbine number
        XCoordinate = zeros(Float64, length(X_vec) , 1, WindFarm.N);               #Array for X coordinates of all points
        YCoordinate = zeros(Float64, length(X_vec) , Real_Rotor_Res, WindFarm.N);  #Array for Y Coordinates of all points
        ZCoordinate = zeros(Float64, 1 , Real_Rotor_Res, 1);                       #Vector containing all height coordinates

        for i in 1:WindFarm.N
        # Create coordinate array of the structure: 1.Dim: Relative turbine, 2.Dim: RotorPoints, 3.Dim: Absolute turbine
            XCoordinate[:, :, i] .= XCoordinate[:, :, i] .+ X_vec .- WindFarm.x_vec[i];
            for j in 1:length(X_vec)
                YCoordinate[j, :, i] .= YCoordinate[j, :, i] .+ Y_vec;
            end
            YCoordinate[:,:,i] .= YCoordinate[:,:,i] .- WindFarm.y_vec[i]
        end
        ZCoordinate[1,:,1] .= Z_vec

    ### Generate grid for convection velocity computation (needed for momentum conserving superposition)
        if WindFarm.Superpos == "Momentum_Conserving"
        # generate basic grid vector 

        if WindFarm.Uc_Res < 4; WindFarm.Uc_Res=4, println("Uc_Res corrected to minimal value 4") end #correction if resolution is chosen too small

        Y_Uc, Z_Uc, WindFarm.Uc_Res = generate_grid_for_Uc(minimum(WindFarm.y_vec)-1, maximum(WindFarm.y_vec)+1, 0.1, 2.0, WindFarm.Uc_Res) #Height normalised -> to be scaled after turbine diameter and hub height read in in "LoadTurbineData"
        
        #Preassign arrays for coordinates (the x coordinates can be reused from the rotor points array)
        Y_for_Uc            =   zeros(Float64, WindFarm.N , WindFarm.Uc_Res, WindFarm.N);
        Z_for_Uc            =   zeros(Float64, 1 , WindFarm.Uc_Res, 1);
        r_for_Uc            =   zeros(WindFarm.N, WindFarm.Uc_Res, WindFarm.N);
        #Preassign arrays for result
        
        
        #Create coordinate arrays
        for ii in 1:WindFarm.N
            # Create coordinate array of the structure: 1.Dim: Relative turbine, 2.Dim: RotorPoints, 3.Dim: Absolute turbine
                for jj in 1:WindFarm.N
                    Y_for_Uc[jj, :, ii] .= Y_Uc .+ WindFarm.y_vec[jj];
                end
                Y_for_Uc[:,:,ii] .= Y_for_Uc[:,:,ii] .- WindFarm.y_vec[ii]
            end
            Z_for_Uc[1,:,1] .= Z_Uc
        
        #Preassignment of additional arrays needed for computation
        u_ambient_for_Uc            = zeros(1, WindFarm.Uc_Res, 1)
        Delta_U_for_Uc              = zeros(WindFarm.N, WindFarm.Uc_Res, WindFarm.N)
        Mixed_wake_for_Uc           = zeros(WindFarm.N, WindFarm.Uc_Res, 1)
        Comutation_Region_ID_for_Uc = trues(WindFarm.N, WindFarm.Uc_Res,WindFarm.N)
        sigma_for_Uc                = zeros(WindFarm.N, WindFarm.Uc_Res, WindFarm.N)
        sigma_m_for_Uc              = zeros(WindFarm.N, WindFarm.Uc_Res, WindFarm.N)
        Lambda_for_Uc               = zeros(1, WindFarm.Uc_Res, WindFarm.N)
        k1_for_Uc                   = zeros(WindFarm.N, WindFarm.Uc_Res, WindFarm.N)  
        k2_for_Uc                   = zeros(WindFarm.N, WindFarm.Uc_Res, WindFarm.N)  
        delta_for_Uc                = zeros(1, WindFarm.Uc_Res, 1)  
        Delta_TI_for_Uc             = zeros(WindFarm.N, WindFarm.Uc_Res, WindFarm.N)
        weighting_Factor_for_Uc     = zeros(WindFarm.N, WindFarm.Uc_Res, WindFarm.N)
        
        else

        #Dummies for struct definition in case other superposition is used 
        Y_for_Uc                    =   zeros(1,1,1)
        Z_for_Uc                    =   zeros(1,1,1) 
        r_for_Uc                    =   zeros(1,1,1)
        u_ambient_for_Uc            =   zeros(1,1,1)
        Delta_U_for_Uc              =   zeros(1,1,1)
        Mixed_wake_for_Uc           =   zeros(1,1,1)
        Comutation_Region_ID_for_Uc =   trues(1,1,1)
        sigma_for_Uc                =   zeros(1,1,1)
        sigma_m_for_Uc              =   zeros(1,1,1)
        Lambda_for_Uc               =   zeros(1,1,1)
        k1_for_Uc                   =   zeros(1,1,1)
        k2_for_Uc                   =   zeros(1,1,1)
        delta_for_Uc                =   zeros(1,1,1)
        Delta_TI_for_Uc             =   zeros(1,1,1)
        weighting_Factor_for_Uc     =   zeros(1,1,1)
        end    
    #= Transform with respect to yaw angle
    TBDone!!
    =#


    # Create struct which holds all computation arrays
    GS=ComputationStruct(WindFarm.name, 0.0, XCoordinate, YCoordinate, ZCoordinate, zeros(length(X_vec) , Real_Rotor_Res, WindFarm.N), Real_Rotor_Res,  alpha_Comp, Yaw_Comp,
                            zeros(1,1,WindFarm.N), zeros(1,1,WindFarm.N), 0, 0, 0, (zeros(1,1,WindFarm.N) .+ WindFarm.u_ambient), (zeros(1,1,WindFarm.N) .+ WindFarm.TI_a),
                            zeros(1,1,WindFarm.N), zeros(1,1,WindFarm.N), zeros(1,1,WindFarm.N), zeros(1,1,WindFarm.N), zeros(1,1,WindFarm.N), zeros(1,1,WindFarm.N), zeros(1,1,WindFarm.N), zeros(1,1,WindFarm.N), 
                            zeros(length(X_vec),Real_Rotor_Res,WindFarm.N), zeros(length(X_vec),Real_Rotor_Res,WindFarm.N), zeros(length(X_vec),Real_Rotor_Res,WindFarm.N), zeros(length(X_vec),Real_Rotor_Res,WindFarm.N), 
                            zeros(1,Real_Rotor_Res,1), zeros(length(X_vec),Real_Rotor_Res,WindFarm.N), trues(length(X_vec),Real_Rotor_Res,WindFarm.N), zeros(length(X_vec),1,WindFarm.N), zeros(length(X_vec),1,1), 
                            zeros(length(X_vec),1,1), zeros(length(X_vec),Real_Rotor_Res,WindFarm.N), Delta_U_for_Uc, Mixed_wake_for_Uc, Y_for_Uc, Z_for_Uc, r_for_Uc, u_ambient_for_Uc, Comutation_Region_ID_for_Uc, 
                            sigma_for_Uc, sigma_m_for_Uc, Lambda_for_Uc, k1_for_Uc, k2_for_Uc, delta_for_Uc, Delta_TI_for_Uc, weighting_Factor_for_Uc, zeros(1,1,WindFarm.N), zeros(1,Real_Rotor_Res,WindFarm.N), 
                            zeros(length(X_vec),Real_Rotor_Res,WindFarm.N), zeros(length(X_vec),Real_Rotor_Res,1), zeros(length(X_vec),Real_Rotor_Res,1), trues(WindFarm.N), 
                            trues(WindFarm.N), falses(WindFarm.N),0, zeros(WindFarm.N), zeros(WindFarm.N), zeros(length(X_vec),1,WindFarm.N)
                        )

    return GS
end#initGraphicArrays

function generate_graphic_grid(WindFarm)
# Generates gridded Array of Y and Z points for graphical computation of the "full flow field" 
    #Create ranges of Y and Z coordinates according to user Input
    y = range(minimum(WindFarm.y_vec)-3, maximum(WindFarm.y_vec)+3, step=WindFarm.Resolution)
    z = range(0, maximum(2 + WindFarm.H/WindFarm.D), length=length(y))
    X = range(minimum(WindFarm.x_vec)-3, maximum(WindFarm.x_vec)+10, step=WindFarm.Resolution)

    #Create gridded matrices of y and z 
    Z, Y = [xi for xi in z, yi in y], [yi for xi in z, yi in y]
    
    # Flatten Array
    Z=vec(Z)
    Y=vec(Y)

    # Find resolution 
    Real_Rotor_Res = length(Y)

    return X, Y, Z, Real_Rotor_Res

end

function correctGS(GS, CS, WindFarm)
# Assigns all already computed values (stored in CS) to graphic struct GS

    # Turbine Data
        GS.Ct_vec = CS.Ct_vec #Ct of each turbine
        GS.P_vec  = CS.P_vec  #P of each turbine
    
    # Coordinates
        # cancelling normalization
        GS.XCoordinates         .= GS.XCoordinates .* WindFarm.D
        GS.YCoordinates         .= GS.YCoordinates .* WindFarm.D
        GS.ZCoordinates         .= GS.ZCoordinates .* WindFarm.D
        GS.r                    .= sqrt.(GS.YCoordinates.^2 .+ (GS.ZCoordinates.-WindFarm.H).^2) # Compute vector in radial & height direction for computation
        # Points for convection velocity
        if WindFarm.Superpos == "Momentum_Conserving"
            GS.Y_for_Uc         .= GS.Y_for_Uc .* WindFarm.D
            GS.Z_for_Uc         .= GS.Z_for_Uc .* WindFarm.D;
            GS.r_for_Uc         .= sqrt.(GS.Y_for_Uc.^2 .+ (GS.Z_for_Uc .- WindFarm.H).^2) # Compute vector in radial & height direction for computation
        end

        GS.i = 1 #Fix iteration counter to 1 (no iteration in graphic comp.)



    return GS

end#correctGS

mutable struct ComputationStruct
#Definition of struct (preassignment of arrays & space)
    name::String;      #Name of the InputFile/ WindFarm used
    CompTime::Float64; #Total computation time (updated after computation)
    #Coordinates/ Arrays
    XCoordinates::Array{Float64,3};
    YCoordinates::Array{Float64,3};
    ZCoordinates::Array{Float64,3};
    RotorPointWeights::Array{Float64,3};
    r::Array{Float64,3}; #Needed for single wake computation. Vector in radial & height direction.
    Real_Rotor_Res::Int; #Amount of computation points per rotor
    # Ambient data
    alpha_Comp::Float64;
    # Turbine specifics
    Yaw_Comp::Vector{Float64};      #Yawangle of each turbine
    Ct_vec::Array{Float64,3};       #Ct of each turbine
    P_vec::Array{Float64,3};        #P of each turbine
    TotalPower::Float64;            #Wind farm total power
    Interp_Ct::Any; #Interpolation function of thrustcoefficient
    Interp_P::Any;  #Interpolation function of power
    u_0_vec::Array{Float64,3};      #Inflow velocity of each turbine (Hubheight)
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
    Computation_Region_ID::BitArray{3}; #ID for limiting computation of the wake region
    #Arrays for Superposition & Meandering
    u_c_vec::Array{Float64,3};          #local convection velocity (for each turbine)
    U_c_Farm::Array{Float64,3};         #Global convection velocity (on farm scale)
    U_c_Farm_old::Array{Float64,3};     #Old global convection velocity from last iteration
    weighting_Factor::Array{Float64,3}; #weighting factor (stores uc_i/U_c_Farm with corrections/ filters)
    #Arrays for computation of U_c_Farm
    Delta_U_for_Uc::Array{Float64,3};        #Velocity deficit needed for computation of global convection velocity 
    Mixed_wake_for_Uc::Array{Float64,3};    #Mixed wake (farm level) for computation of global convection velocity
    Y_for_Uc::Array{Float64,3};
    Z_for_Uc::Array{Float64,3};
    r_for_Uc::Array{Float64,3};
    u_ambient_for_Uc::Array{Float64,3};
    Computation_Region_ID_for_Uc::BitArray{3};#ID for limiting computation of the wake region
    sigma_for_Uc::Array{Float64,3};
    sigma_m_for_Uc::Array{Float64,3};
    Lambda_for_Uc::Array{Float64,3};
    k1_for_Uc::Array{Float64,3};
    k2_for_Uc::Array{Float64,3};
    delta_for_Uc::Array{Float64,3};
    Delta_TI_for_Uc::Array{Float64,3};
    weighting_Factor_for_Uc::Array{Float64,3};
    #Arrays exclusively for Meandering
    psi::Array{Float64,3};      #fluctuation intensitys
    Lambda::Array{Float64,3};   #Integral length scale of representative eddy
    sigma_m::Array{Float64,3};  #wake width correction parameter
    #Arrays needed to superimpose
    U_Farm::Array{Float64,3};       #Final superimposed result for Velocity
    TI_Farm::Array{Float64,3};      #Final superimposed result for turbulence intensity
    #Computation Parameters
    ID_OutOperConst::BitVector;     # Identification of turbines which are out of operation
    ID_Turbines::BitVector;         # Identification of turbines which should not be computed in following iteration (already converged)
    ID_Turbines_Computed::BitVector;# Turbines whose wake have already been computed
    i::Int;                         # iteration counter
    StreamwiseOrder::Vector{Int64}; # Vector to identify streamwise order of turbines
    CompOrder::Vector{Int64};       # Vector to identify computation order/ priority
    #Temporary computation help
    tmp::Array{Float64,3};
end #mutable struct "ComputationStruct"

mutable struct GraphicComputationStruct
    #Definition of struct (preassignment of arrays & space)
        name::String;      #Name of the InputFile/ WindFarm used
        CompTime::Float64; #Total computation time (updated after computation)
        #Coordinates/ Arrays
        XCoordinates::Array{Float64,3};
        YCoordinates::Array{Float64,3};
        ZCoordinates::Array{Float64,3};
        RotorPointWeights::Array{Float64,3};
        r::Array{Float64,3}; #Needed for single wake computation. Vector in radial & height direction.
        Real_Rotor_Res::Int; #Amount of computation points per rotor
        # Ambient data
        alpha_Comp::Float64;
        # Turbine specifics
        Yaw_Comp::Vector{Float64};      #Yawangle of each turbine
        Ct_vec::Array{Float64,3};       #Ct of each turbine
        P_vec::Array{Float64,3};        #P of each turbine
        TotalPower::Float64;            #Wind farm total power
        Interp_Ct::Any; #Interpolation function of thrustcoefficient
        Interp_P::Any;  #Interpolation function of power
        u_0_vec::Array{Float64,3};      #Inflow velocity of each turbine (Hubheight)
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
        Computation_Region_ID::BitArray{3}; #ID for limiting computation of the wake region
        #Arrays for Superposition & Meandering
        u_c_vec::Array{Float64,3};          #local convection velocity (for each turbine)
        U_c_Farm::Array{Float64,3};         #Global convection velocity (on farm scale)
        U_c_Farm_old::Array{Float64,3};     #Old global convection velocity from last iteration
        weighting_Factor::Array{Float64,3}; #weighting factor (stores uc_i/U_c_Farm with corrections/ filters)
        #Arrays for computation of U_c_Farm
        Delta_U_for_Uc::Array{Float64,3};        #Velocity deficit needed for computation of global convection velocity 
        Mixed_wake_for_Uc::Array{Float64,3};    #Mixed wake (farm level) for computation of global convection velocity
        Y_for_Uc::Array{Float64,3};
        Z_for_Uc::Array{Float64,3};
        r_for_Uc::Array{Float64,3};
        u_ambient_for_Uc::Array{Float64,3};
        Computation_Region_ID_for_Uc::BitArray{3};#ID for limiting computation of the wake region
        sigma_for_Uc::Array{Float64,3};
        sigma_m_for_Uc::Array{Float64,3};
        Lambda_for_Uc::Array{Float64,3};
        k1_for_Uc::Array{Float64,3};
        k2_for_Uc::Array{Float64,3};
        delta_for_Uc::Array{Float64,3};
        Delta_TI_for_Uc::Array{Float64,3};
        weighting_Factor_for_Uc::Array{Float64,3};
        #Arrays exclusively for Meandering
        psi::Array{Float64,3};      #fluctuation intensitys
        Lambda::Array{Float64,3};   #Integral length scale of representative eddy
        sigma_m::Array{Float64,3};  #wake width correction parameter
        #Arrays needed to superimpose
        U_Farm::Array{Float64,3};       #Final superimposed result for Velocity
        TI_Farm::Array{Float64,3};      #Final superimposed result for turbulence intensity
        #Computation Parameters
        ID_OutOperConst::BitVector;     # Identification of turbines which are out of operation
        ID_Turbines::BitVector;         # Identification of turbines which should not be computed in following iteration (already converged)
        ID_Turbines_Computed::BitVector;# Turbines whose wake have already been computed
        i::Int;                         # iteration counter
        StreamwiseOrder::Vector{Int64}; # Vector to identify streamwise order of turbines
        CompOrder::Vector{Int64};       # Vector to identify computation order/ priority
        #Temporary computation help
        tmp::Array{Float64,3};
end #mutable struct "ComputationStruct"

end #Module