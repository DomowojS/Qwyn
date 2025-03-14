#=
Postprocessing module
- Compute Graphical output data
- Plot as requested
=#
module Postprocessing
using PlotlyJS: Plot, scatter, Layout, attr
using GLMakie: Figure, Axis, contourf!, Colorbar

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

function AdvancedPlots(WindFarm, GS, Levels)

    #Create grids for Plotting
    x = GS.XCoordinates[:, 1, 1]./WindFarm.D
    y = GS.YCoordinates[1, :, 1]./WindFarm.D
    x_grid = repeat(x, 1, length(y))
    y_grid = repeat(y', length(x), 1)


    # Rotate system, so wind farm always looks identical (270° inflow left to right)
    rotation_angle = deg2rad(270 - WindFarm.alpha)
    x_grid_rotated = x_grid .* cos(rotation_angle) .- y_grid .* sin(rotation_angle)
    y_grid_rotated = x_grid .* sin(rotation_angle) .+ y_grid .* cos(rotation_angle)


    # Read wind data
    z = GS.U_Farm[:, :, 1]./WindFarm.u_ambient_zprofile_4Graphic[1,1,1]

    # Create a figure
    fig = Figure(size = (800, 600))
    
    # Create an axis
    ax = Axis(
        fig[1, 1],
        xlabel = "x/D",
        ylabel = "y/D",
    )
    
    # Create the contour plot
    co=contourf!(ax, x_grid_rotated, y_grid_rotated, z, colormap = :inferno, levels=range(0, 1, length = Levels))
    
    # Add a colorbar
    Colorbar(fig[1, 2], co)
    
    # Display the figure
    display(fig)
    
end

end#Postprocessing Module