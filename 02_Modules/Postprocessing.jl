#=
Postprocessing module
- Safe results in result struct
- Compute Graphical output data
- Plot as requested
=#
module Postprocessing
using PlotlyJS
export SimplePlots, AdvancedPlots

function SimplePlots(WindFarm, CS)
#This function plots all simple requested plots from the input file    
    if WindFarm.Plot_power == true
    # Power output plot
        # Extract the indication
        indication = WindFarm.Turbine_Identification 
        # Create an evenly spaced vector for x-axis
        x = 1:length(indication)

        # Get the corresponding entries from CS.P_vec
        if WindFarm.Normalize_to == 0
            y = CS.P_vec[indication]./1000000
            min, = findmin(y)
            max, = findmax(y)
            dif= max - min
            min = min- 0.05dif
            max = max+ 0.05dif
                # Create the plot
            plot = Plot(
                scatter(
                    x = x,
                    y = y,
                    mode = "lines+markers",
                    name = "P_vec"
                ),
                Layout(
                    title = "Modelled Power Output",
                    xaxis = attr(title = "Turbine number"),
                    yaxis = attr(title = "P in MW"),
                    range = [min, max]
                )
            )
        else
            y = CS.P_vec[indication]./CS.P_vec[WindFarm.Normalize_to]
            min, = findmin(y)
            max, = findmax(y)
            dif= max - min
            #min = min- 0.05dif
            #max = max+ 0.05dif
            min = 0
            max = 1
                # Create the plot
                plot = Plot(
                    scatter(
                        x = x,
                        y = y,
                        mode = "lines+markers",
                        name = "P_vec"
                    ),
                    Layout(
                        title = "Modelled Power Output",
                        xaxis = attr(
                        title = "Turbine number",
                        tickmode = "array",
                        tickvals = x,
                        ticktext = string.(x),
                        titlefont = attr(size = 18, family = "Arial, sans-serif")
                        ),
                        yaxis = attr(
                            title = "P/P₁",
                            titlefont = attr(size = 18, family = "Arial, sans-serif"),
                            range = [min, max]
                        )
                    )
                )
        end
        # Display the plot
    display(plot)
    end
    if WindFarm.Plot_windspeed == true
    # Wind speed plot
        # Extract the indication
        indication = WindFarm.Turbine_Identification 
        # Create an evenly spaced vector for x-axis
        x = 1:length(indication)

        # Get the corresponding entries from CS.u_0_vec
        if WindFarm.Normalize_to == 0
            y = CS.u_0_vec[indication]
            min, = findmin(y)
            max, = findmax(y)
            dif= max - min
            min = min- 0.05dif
            max = max+ 0.05dif
                # Create the plot
            plot1 = Plot(
                scatter(
                    x = x,
                    y = y,
                    mode = "lines+markers",
                    name = "u_0_vec"
                ),
                Layout(
                    title = "Modelled wind speed at hub height",
                    xaxis = attr(title = "Turbine number"),
                    yaxis = attr(title = "u₀ in m/s"),
                    range = [min, max]
                )
            )
        else
            y = CS.u_0_vec[indication]./CS.u_0_vec[WindFarm.Normalize_to]
            min, = findmin(y)
            max, = findmax(y)
            dif= max - min
            min = min- 0.05dif
            max = max+ 0.05dif
                # Create the plot
                plot1 = Plot(
                    scatter(
                        x = x,
                        y = y,
                        mode = "lines+markers",
                        name = "u_0_vec"
                    ),
                    Layout(
                        title = "Modelled wind speed at hub height",
                        xaxis = attr(
                        title = "Turbine number",
                        tickmode = "array",
                        tickvals = x,
                        ticktext = string.(x),
                        titlefont = attr(size = 18, family = "Arial, sans-serif")
                        ),
                        yaxis = attr(
                            title = "uᵢ/u₁",
                            titlefont = attr(size = 18, family = "Arial, sans-serif"),
                            range = [min, max]
                        )
                    )
                )
        end
        # Display the plot
    display(plot1)
    end
    if WindFarm.Plot_turbulence == true
    # Turbulence intensity plot
        # Extract the indication
        indication = WindFarm.Turbine_Identification 
        # Create an evenly spaced vector for x-axis
        x = 1:length(indication)

        # Get the corresponding entries from CS.TI_0_vec
        if WindFarm.Normalize_to == 0
            y = CS.TI_0_vec[indication]
            min, = findmin(y)
            max, = findmax(y)
            dif= max - min
            min = min- 0.05dif
            max = max+ 0.05dif
                # Create the plot
            plot2 = Plot(
                scatter(
                    x = x,
                    y = y,
                    mode = "lines+markers",
                    name = "TI_0_vec"
                ),
                Layout(
                    title = "Modelled turbulence intensity at hub height",
                    xaxis = attr(title = "Turbine number"),
                    yaxis = attr(title = "TI₀ in %"),
                    range = [min, max]
                )
            )
        else
            y = CS.TI_0_vec[indication]./CS.TI_0_vec[WindFarm.Normalize_to]
            min, = findmin(y)
            max, = findmax(y)
            dif= max - min
            min = min- 0.05dif
            max = max+ 0.05dif
                # Create the plot
                plot2 = Plot(
                    scatter(
                        x = x,
                        y = y,
                        mode = "lines+markers",
                        name = "TI_0_vec"
                    ),
                    Layout(
                        title = "Modelled turbulence Intensity at hub height",
                        xaxis = attr(
                        title = "Turbine number",
                        tickmode = "array",
                        tickvals = x,
                        ticktext = string.(x),
                        titlefont = attr(size = 18, family = "Arial, sans-serif")
                        ),
                        yaxis = attr(
                            title = "TIᵢ/TI₁",
                            titlefont = attr(size = 18, family = "Arial, sans-serif"),
                            range = [min, max]
                        )
                    )
                )
        end
        # Display the plot
    display(plot2)
    end

end#SimplePlots

function AdvancedPlots(WindFarm, GS)
#This function plots the advances plots (flow fields) 
    
    if WindFarm.Plot_wind_field == true
        z_level = WindFarm.H  # Replace with your specific Z-coordinate value

        # Find the closest z-level in GS.ZCoordinates
        z_coords = GS.ZCoordinates[1, :, 1]  # Extract relevant z values from GS.ZCoordinates
        z_level = z_coords[argmin(abs.(z_coords .- z_level))]  # Closest z-level

        # Extract the relevant slices
        x = GS.XCoordinates[:, 1, 1]                  # X coordinates (size 55)
        y = GS.YCoordinates[1, :, 1]                  # Y coordinates (size 729)
        u_values = GS.U_Farm[:, :, 1]                 # U_Farm values for contour plot

        # Filter data based on the specified z level
        y_filtered = y[z_coords .== z_level]            # Only y values where z matches z_level
        u_filtered = u_values[:, z_coords .== z_level]  # Corresponding U_Farm values

        # Create the 2D contour plot
        plt = plot(
            contour(
                y = x./WindFarm.D,
                x = y_filtered./WindFarm.D,
                z = transpose(u_filtered),
                colorscale = "Viridis",
                contours_coloring = "heatmap",
                colorbar_title = "U_Farm",
                linewidth = 0  # Set contour line width to 0
            ),
            Layout(
                xaxis_title = "x/D",  # X-axis label
                yaxis_title = "y/D"   # Y-axis label
            )
        )
    
        # Display the plot
        display(plt)
    end

end#AdvancedPlots

end#Postprocessing Module