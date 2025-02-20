#include("02_Modules/PKG_Manager.jl")            #pck manager (for own overiew -> !!to be deleted!! 
include("02_Modules/Input_Processing.jl")       #Module to process input data
include("02_Modules/Initialisation_Module.jl")  #Module for array initialisation/ space preallocation
include("02_Modules/SimpleComputation.jl")      #Module for simple computation.
include("02_Modules/Optimisation.jl")      #Module for optimisation
include("02_Modules/Postprocessing.jl")
using .Input_Processing, .Initialisation_Module, .SimpleComputation, .Optimisation, .Postprocessing, MAT, Base.Threads, PlotlyJS, JLD2, FileIO


function Qwyn_Simple(u_ambient::Real, alpha::Real, TI_a::Real)
#=  This function excexutes Qwyn. 
    Input configs are taken from 01_Inputs.
    All function blocks can be found in 02_Modules
=#

    t_start = time(); #Timer
    WF = generateWF("WF", "01_Inputs", u_ambient, alpha, TI_a) #Generate array consisting of all input data for each Input 
                                                               #file in "01_Inpts"

    Results=Vector(undef, length(WF)); #Preassign result struct
    
    # Iterate over defined cases. Compute one by one (according to settings)
    i=0;#Loop counte
    for WindFarm in WF 
    t_start_loop = time();#Looptimer       
    i=i+1;
        println("###########################")
        println("Computing: ", WindFarm.name)   #Terminal output for which input file is being processed
        println("###########################")

        #Initialise all arrays & matrices needed for the computation.
        println("Initialising arrays..")
        WindFarm, CS = initCompArrays(WindFarm) #Initialises mutable struct "CA" which contains necessary 
                                                #computation arrays & computes coordinates acc. to user Input.

        println("..loading turbine data..")
        LoadTurbineDATA!(WindFarm, CS)          #Update Input & computation structs with provided power & 
                                                #thrust curves
        
        println("..loading atmospheric data..")
        WindFarm.u_ambient_zprofile, CS = LoadAtmosphericData(WindFarm,CS)      #Update Input & computation structs with atmospheric data &
                                                                                #(wind shear profile, wind rose etc.)
        
        println("..finding computation order..")
        FindStreamwiseOrder!(WindFarm, CS)



        # Iterating over turbine rows
        println("Starting computation...")
        for CS.i=1:maximum(CS.CompOrder)
            if WindFarm.Superpos=="Momentum_Conserving" println("Step ", CS.i) end

            FindComputationRegion!(WindFarm, CS)
            
            if any(CS.ID_Turbines != 0) # Check if any turbine's wake has to be computed

            Single_Wake_Computation!(WindFarm, CS)  #Compute single wake effect
            
            Superposition!(WindFarm, CS, WindFarm.u_ambient_zprofile)            #Compute mixed wake
            
            getTurbineInflow!(WindFarm, CS)         #Evaluate new inflow data
            
            getNewThrustandPower!(WindFarm, CS)     #Evaluate new operation properties

            end
        end
            
        getTotalPower!(CS)  #Compute total power of the wind farm 
        println("...computation finished!")
        t_end_loop = time()-t_start_loop;
        println("Computation time of input No.$i: $(round(t_end_loop, digits = 5)) s.")
        for i in 1:2 println(".") end

        if any([WindFarm.Plot_power, WindFarm.Plot_windspeed, WindFarm.Plot_turbulence, WindFarm.Plot_wind_field, WindFarm.Plot_turbulence_field]) == true
        println("Plots on the way...")
            if any([WindFarm.Plot_power, WindFarm.Plot_windspeed, WindFarm.Plot_turbulence]) == true
                SimplePlots(WindFarm, CS)
            end

            if any([WindFarm.Plot_wind_field, WindFarm.Plot_turbulence_field]) == true
                GS = Qwyn_Graphic(WindFarm, CS)
                AdvancedPlots(WindFarm, GS)
            end

            println("...finished!")
        end
    
     #TMP=reshape(CS.P_vec[[4, 12, 20, 28, 36, 44, 52, 60]]./CS.P_vec[4], 8) #270
     #TMP=reshape(CS.P_vec[[5, 12, 19, 26, 33]]./CS.P_vec[4], 5)              #222
     #TMP=reshape(CS.P_vec[[4, 13, 22, 31, 40]]./CS.P_vec[4], 5)             #312
     #Base.print_matrix(stdout, TMP)
        
        # Organize final output
        if WindFarm.Extended_Output == true #If extended output is requested 
            CS.CompTime=t_end_loop;
            Results[i]=CS;
        elseif WindFarm.Extended_Output == false
            Results[i]=ShortResult(WindFarm.name, t_end_loop, WindFarm.x_vec, WindFarm.y_vec, vec(WindFarm.u_ambient_zprofile), 
                        vec(CS.Ct_vec), vec(CS.P_vec), CS.TotalPower, vec(CS.u_0_vec), vec(CS.TI_0_vec), CS.U_Farm, CS.TI_Farm);
        end
    end#Loop over Input structs

        println("Total processing time of $i case(s): $(round(time()-t_start, digits = 4))s.")
        return Results, WF #Returns results & Inputfiles as structs

end#Qwyn_Simple

###############################################################################
# TEMPORARY ONLY FOR STUDY #

function Qwyn_AEP_RPStudy(TI_a::Real ,path2windrose, RPoints)
    #=  This functions runs an AEP computation. 
        Input configs are taken from 01_Inputs.
        All function blocks can be found in 02_Modules
    
        Inputs:
        1) Turbulence Intensity
        2) Path to wind rose .txt file
    =#
        t_start = time(); #Timer
    
        #Read provided wind rose data
        angles, speeds, frequencies=read_Windrose_data(path2windrose)
    
        WF = generateWF("WF", "01_Inputs", speeds[1], angles[1], TI_a) #Generate array consisting of all input data for each Input 
                                                                       #file in "01_Inpts"
        
        AEP_Results=Vector(undef, length(WF));     #Preassign AEP-result vector
        Powers_AllCases=Vector{Float64}(undef, length(angles))
        TotalAEP=rand(1,1)
        TotalAEP_SingleTurbine=Vector{Float64}(undef, 80)

        # Iterate over defined cases. Compute one by one (according to settings)
        i=0;#Loop counter
        for WindFarm in WF 
        t_start_loop = time();#Looptimer       
        i=i+1;
            println("###########################")
            println("Number of Rotor Points: $RPoints")   #Terminal output for which input file is being processed
            println("###########################")
    
            
        P_vec_AllCases=Vector{Array{Float64,3}}(); #Reset P_vec from Config to Config
    
    
        ######## START LOOP OVER all wind speeds HERE ############################
            
            last_printed_percent = Threads.Atomic{Int}(0)  # Atomic counter for progress
            WindFarm.Rotor_Res = RPoints # Rotor Resolution
            
            for case in 1:length(frequencies)
    
                # Assign right wind speed & angle for this iteration
                WindFarm.u_ambient = speeds[case] #speed
                WindFarm.alpha     = angles[case] #angle
                
                #Initialise all arrays & matrices needed for the computation.
                WindFarm, CS = initCompArrays(WindFarm) #Initialises mutable struct "CA" which contains necessary 
                                                        #computation arrays & computes coordinates acc. to user Input.
    
                LoadTurbineDATA!(WindFarm, CS)          #Update Input & computation structs with provided power & 
                                                        #thrust curves
                
                WindFarm.u_ambient_zprofile, CS = LoadAtmosphericData(WindFarm,CS)      #Update Input & computation structs with atmospheric data &
                                                                                        #(wind shear profile, wind rose etc.)
                
                FindStreamwiseOrder!(WindFarm, CS)
    
    
    
                # Iterating over turbine rows
                for CS.i=1:maximum(CS.CompOrder)
                    FindComputationRegion!(WindFarm, CS)
                    
                    if any(CS.ID_Turbines != 0) # Check if any turbine's wake has to be computed
    
                    Single_Wake_Computation!(WindFarm, CS)  #Compute single wake effect
                    
                    Superposition!(WindFarm, CS, WindFarm.u_ambient_zprofile)            #Compute mixed wake
                    
                    getTurbineInflow!(WindFarm, CS)         #Evaluate new inflow data
                    
                    getNewThrustandPower!(WindFarm, CS)     #Evaluate new operation properties
    
                    end
                end
                    
                getTotalPower!(CS)  #Compute total power of the wind farm 
    
                Powers_AllCases[case]=CS.TotalPower;   #Save computed power for this case
                push!(P_vec_AllCases, CS.P_vec)           #Save computed power for each single turbine for this case
    
                 # Calculate and display progress every full percentage
                    progress = round(Int, 100 * case / length(frequencies))
                    old_progress = Threads.atomic_cas!(last_printed_percent, progress-1, progress)
                    if old_progress != progress
                        println("Progress: $progress%")
                    end
    
            end
                
                TotalAEP=sum(Powers_AllCases.*frequencies) * 8760;               #Compute total AEP by multyplying power with frequencies
                TotalAEP_SingleTurbine=sum(P_vec_AllCases.*frequencies) * 8760; #Compute total AEP by multyplying power with frequencies
    
               # Display end of computation 
                println("Computation finished")
                println("###########################")
        end#Loop over Input structs
    
       return TotalAEP, TotalAEP_SingleTurbine
       #return TotalAEP
    
    end#Qwyn_AEP_RPStudy

###############################################################################

function Qwyn_AEP(TI_a::Real ,path2windrose)
#=  This functions runs an AEP computation. 
    Input configs are taken from 01_Inputs.
    All function blocks can be found in 02_Modules

    Inputs:
    1) Turbulence Intensity
    2) Path to wind rose .txt file
=#
    t_start = time(); #Timer

    #Read provided wind rose data
    angles, speeds, frequencies=read_Windrose_data(path2windrose)

    WF = generateWF("WF", "01_Inputs", speeds[1], angles[1], TI_a) #Generate array consisting of all input data for each Input 
                                                                   #file in "01_Inpts"

    AEP_Results=Vector(undef, length(WF));     #Preassign AEP-result vector
    Powers_AllCases=Vector{Float64}(undef, length(angles))

    # Iterate over defined cases. Compute one by one (according to settings)
    i=0;#Loop counter
    for WindFarm in WF 
    t_start_loop = time();#Looptimer       
    i=i+1;
        println("###########################")
        println("Computing: ", WindFarm.name)   #Terminal output for which input file is being processed
        println("###########################")

        
    P_vec_AllCases=Vector{Array{Float64,3}}(); #Reset P_vec from Config to Config


    ######## START LOOP OVER all wind speeds HERE ############################
        
        last_printed_percent = Threads.Atomic{Int}(0)  # Atomic counter for progress
        
        for case in 1:length(frequencies)

            # Assign right wind speed & angle for this iteration
            WindFarm.u_ambient = speeds[case] #speed
            WindFarm.alpha     = angles[case] #angle

            #Initialise all arrays & matrices needed for the computation.
            WindFarm, CS = initCompArrays(WindFarm) #Initialises mutable struct "CA" which contains necessary 
                                                    #computation arrays & computes coordinates acc. to user Input.

            LoadTurbineDATA!(WindFarm, CS)          #Update Input & computation structs with provided power & 
                                                    #thrust curves
            
            WindFarm.u_ambient_zprofile, CS = LoadAtmosphericData(WindFarm,CS)      #Update Input & computation structs with atmospheric data &
                                                                                    #(wind shear profile, wind rose etc.)
            
            FindStreamwiseOrder!(WindFarm, CS)



            # Iterating over turbine rows
            for CS.i=1:maximum(CS.CompOrder)
                FindComputationRegion!(WindFarm, CS)
                
                if any(CS.ID_Turbines != 0) # Check if any turbine's wake has to be computed

                Single_Wake_Computation!(WindFarm, CS)  #Compute single wake effect
                
                Superposition!(WindFarm, CS, WindFarm.u_ambient_zprofile)            #Compute mixed wake
                
                getTurbineInflow!(WindFarm, CS)         #Evaluate new inflow data
                
                getNewThrustandPower!(WindFarm, CS)     #Evaluate new operation properties

                end
            end
                
            getTotalPower!(CS)  #Compute total power of the wind farm 

            Powers_AllCases[case]=CS.TotalPower;   #Save computed power for this case
            push!(P_vec_AllCases, CS.P_vec)           #Save computed power for each single turbine for this case

             # Calculate and display progress every full percentage
                progress = round(Int, 100 * case / length(frequencies))
                old_progress = Threads.atomic_cas!(last_printed_percent, progress-1, progress)
                if old_progress != progress
                    println("Progress: $progress%")
                end

        end
            
            TotalAEP=sum(Powers_AllCases.*frequencies) * 8760;               #Compute total AEP by multyplying power with frequencies
            TotalAEP_SingleTurbine=sum(P_vec_AllCases.*frequencies) * 8760; #Compute total AEP by multyplying power with frequencies
        
            t_end_loop = time()-t_start_loop;


            AEP_Results[i] = AEP(WindFarm.name, t_end_loop, TotalAEP, Powers_AllCases,  reshape(TotalAEP_SingleTurbine, WindFarm.N), P_vec_AllCases, angles, speeds, frequencies);

           # Display end of computation 
            println("Computation finished")
            println("Computation time of input No.$i: $(round(t_end_loop, digits = 5)) s.")
            println("###########################")
    end#Loop over Input structs

    println("Total processing time of $i case(s): $(round(time()-t_start, digits = 4))s.")
    return AEP_Results #Returns results & Inputfiles as structs

end#Qwyn_AEP

function Qwyn_Layout_Optimiser(TI_a::Real ,path2windrose)
#=  This functions optimises the turbine coordinates for higher AEP. 
    Inputs configs are taken from 01_Inputs.
    All function blocks can be found in 02_Modules

    Inputs:
    1) Turbulence Intensity
    2) Path to wind rose .txt file
=#

    t_start = time(); #Timer

    #Read provided wind rose data
    angles, speeds, frequencies=read_Windrose_data(path2windrose)

    WF = generateWF("WF", "01_Inputs", speeds[1], angles[1], TI_a) #Generate array consisting of all input data for each Input 
                                                                   #file in "01_Inpts"

    #Preassign result vectors
    Layout_Optim_Results=Vector(undef, length(WF));                #Preassign Optim-result vector
    Powers_AllCases=Vector{Float64}(undef, length(angles))         #Preassign Temporary result vector used within Optim target fct

    # Iterate over defined cases. Compute one by one (according to settings)
    i=0;#Loop counter
    for WindFarm in WF 
    t_start_loop = time();#Looptimer       
    i=i+1;
        println("###########################")
        println("Optimising: ", WindFarm.name)   #Terminal output for which input file is being processed
        println("###########################")
        OptimalAEP, x_vec_optimal, y_vec_optimal, ret, num_evals = SetAndRunOptimiser(WindFarm, speeds, angles, frequencies, Powers_AllCases);
        
        # Give some useroutpuy
        println(
            """
            Optimal AEP             : $OptimalAEP
            Optimal x coordinates   : $x_vec_optimal
            Optimal y coordinates   : $y_vec_optimal
            solution status         : $ret
            # function evaluation   : $num_evals
            """
        )

        # Quick & dirty machine position plot
        p = plot_turbine_positions(x_vec_optimal, y_vec_optimal)
        display(p)

        # Update timer
        t_end_loop = time()-t_start_loop;

        # Display end of computation 
        println("Optimisation finished")
        println("Optimisation time of input No.$i: $(round(t_end_loop, digits = 5)) s.")
        println("###########################")
            
        # Write to result struct
        Layout_Optim_Results[i] = Layout_Optim_Results_Struct(WindFarm.name, t_end_loop, OptimalAEP, WindFarm.x_vec, WindFarm.y_vec, x_vec_optimal, y_vec_optimal);

    ########## Add normal "AEP" function call here --> For using 


    end#Loop over Input structs


    println("Total processing time of $i case(s): $(round(time()-t_start, digits = 4))s.")
    return Layout_Optim_Results #Returns results & Inputfiles as structs                                                       

end#Qwyn_Optimiser

function Qwyn_Graphic(WindFarm, CS)
# Computes full flowfield at resolution as specified in InputFile
    println("...computing Full Flow Field...")

    GS = initGraphicArrays(WindFarm)        #Initialises mutable struct "GS" which contains necessary 
                                            #computation arrays for Graphical Output
    
    # Correct preassigned values in "GS" according to already computed data in "CS"
    GS = correctGS(GS, CS, WindFarm)

                                  
    WindFarm.u_ambient_zprofile_4Graphic, GS = LoadAtmosphericData(WindFarm, GS)       #Update Input & computation structs with atmospheric data &
                                                                                       #(wind shear profile, wind rose etc.)                                       
    
    # Set Comp & Streamwise order to 1. No iteration needed, only one parallel computation for all turbines
    GS.CompOrder       .= 1
    GS.StreamwiseOrder .= 1

    if WindFarm.Superpos=="Momentum_Conserving" println("Step ", CS.i) end

        FindComputationRegion!(WindFarm, GS)

        Single_Wake_Computation!(WindFarm, GS)  #Compute single wake effect
    
        Superposition!(WindFarm, GS, WindFarm.u_ambient_zprofile_4Graphic)            #Compute mixed wake

    return GS #Returns Graphical Struct

end#Qwyn_Graphic

function plot_turbine_positions(x_vec_optimal, y_vec_optimal)
# Plots turbine positions with boundaries after optimisation

    # Create hover text for each point with correct HTML formatting
    hover_text = [string("Turbine ", i, "<br>", 
                        "x: ", round(x_vec_optimal[i], digits=2), "<br>",
                        "y: ", round(y_vec_optimal[i], digits=2))
                 for i in 1:length(x_vec_optimal)]
    
    # Create the scatter trace for turbines
    turbine_trace = scatter(
        x = x_vec_optimal,
        y = y_vec_optimal,
        mode = "markers+text",
        marker = attr(
            size = 15,
            color = "rgb(31, 119, 180)",
            line = attr(color = "rgb(0, 0, 0)", width = 1)
        ),
        text = ["T$i" for i in 1:length(x_vec_optimal)],
        textposition = "top center",
        hovertext = hover_text,
        hoverinfo = "text",
        name = "Turbines"
    )
    
    # Create boundary box
    boundary_trace = scatter(
        #x = [minimum(x_vec_optimal), maximum(x_vec_optimal), maximum(x_vec_optimal), minimum(x_vec_optimal), minimum(x_vec_optimal)],  # Close the box by returning to start
        #y = [minimum(y_vec_optimal), minimum(y_vec_optimal), maximum(y_vec_optimal), maximum(y_vec_optimal), minimum(y_vec_optimal)],
        #For Horns Rev
        x = [6, 69, 63, 0, 6],  # Close the box by returning to start
        y = [-48.63, -48.63, 0, 0, -48.63],
        mode = "lines",
        line = attr(
            color = "red",
            width = 2,
            dash = "dash"
        ),
        name = "Outer Boundaries"
    )
    
    # Create circles for minimum distance
    traces = GenericTrace{Dict{Symbol, Any}}[turbine_trace, boundary_trace]
    
    for i in 1:length(x_vec_optimal)
        # Generate points for a circle
        t = range(0, 2π, length=100)
        circle_x = x_vec_optimal[i] .+ 3.7 .* cos.(t)
        circle_y = y_vec_optimal[i] .+ 3.7 .* sin.(t)
        
        circle_trace = scatter(
            x = circle_x,
            y = circle_y,
            mode = "lines",
            line = attr(
                color = "rgba(169, 169, 169, 0.3)",  # Semi-transparent gray
                width = 1
            ),
            showlegend = i == 1,  # Only show in legend once
            name = "Min Distance (3.7)",
            hoverinfo = "skip"
        )
        push!(traces, circle_trace)
    end
    
    # Create the layout
    layout = Layout(
        title = "Optimal Turbine Positions",
        xaxis = attr(
            title = "X Coordinate",
            zeroline = true,
            gridcolor = "rgb(235, 235, 235)",
            showgrid = true,
            #ange = [minimum(x_vec_optimal)-1, maximum(x_vec_optimal)+1]  # Slightly larger than boundaries for visibility
            # For HR1
            range = [0-1, 69+1]  # Slightly larger than boundaries for visibility
        ),
        yaxis = attr(
            title = "Y Coordinate",
            zeroline = true,
            gridcolor = "rgb(235, 235, 235)",
            showgrid = true,
            scaleanchor = "x",
            scaleratio = 1,
            #range = [minimum(y_vec_optimal)-1, maximum(y_vec_optimal)+1]  # Slightly larger than boundaries for visibility
            # For HR1
            range = [-48.63-1, 0+1]  # Slightly larger than boundaries for visibility
        ),
        hovermode = "closest",
        showlegend = true,
        legend = attr(
            x = 1.1,
            y = 1.0
        ),
        plot_bgcolor = "rgb(250, 250, 250)",
        width = 900,  # Increased width to accommodate legend
        height = 600
    )
    
    # Create and return the plot with all traces
    p = plot(traces, layout)
    
    return p
end
########## Functional elements ############

# Consice result struct
struct ShortResult
# Short result struct including important inputs
    name::String;
    CompTime::Float64;
    x_Coordinates_streamwise::Vector{Float64};
    y_Coordinates_streamwise::Vector{Float64};
    u_heightprofile::Vector{Float64};

    #Results Turbines:
    Ct_vec::Vector{Float64};
    P_vec::Vector{Float64};
    TotalPower::Float64;
    u_0_av::Vector{Float64};
    TI_0_av::Vector{Float64};
    U_0_RotorDistr::Array{Float64, 3}
    TI_0_RotorDistr::Array{Float64, 3}
end#ShortResult

# Struct for AEP Results including Inputs
struct AEP
# Struct for AEP Resutls
    name::String;
    CompTime::Float64;
    Total_AEP::Float64;
    Powers_AllCases::Vector{Float64};
    TotalAEP_SingleTurbine::Vector{Float64};
    P_vec_AllCases::Vector{Array{Float64, 3}};
    angles::Vector{Float64};
    speeds::Vector{Float64};
    frequencies::Vector{Float64};
end

struct Layout_Optim_Results_Struct
    # Struct for AEP Resutls
        name::String;
        CompTime::Float64;
        Optimised_AEP::Float64;
        x_coordinates_initial::Vector{Float64};
        y_coordinates_initial::Vector{Float64};
        x_coordinates_optimal::Vector{Float64};
        y_coordinates_optimal::Vector{Float64};
end



























                                #= Temporary stuff -> for saving to MAT and plotting in Matlab!    
                                    # Convert the struct to a dictionary
                                    CS.Interp_Ct=0
                                    CS.Interp_P=0
                                struct_dict = Dict{String, Any}(string.(propertynames(CS)) .=> getfield.(Ref(CS), propertynames(CS)))
                                struct_dict2 = Dict{String, Any}(string.(propertynames(WindFarm)) .=> getfield.(Ref(WindFarm), propertynames(WindFarm)))

                                # Specify the filename for the .mat file
                                filename = "99_PlotWMATLAB/WindFarmCS_NewCoordinateSys.mat"
                                filename2 = "99_PlotWMATLAB/WindFarmWF_NewCoordinateSys.mat"
                                # Save the struct to the .mat file
                                matwrite(filename, struct_dict)
                                matwrite(filename2, struct_dict2)
                                =#