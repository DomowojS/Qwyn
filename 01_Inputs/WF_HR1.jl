#=
This is a first Input file
=#
@with_kw mutable struct WF
    # Name of the wind Farm
    name::String="HornsRev1";
    
    # Wind Farm Data
    N::Int=6; #Number of turbines
    
    x_vec=[0, 7, 14, 21, 28, 35]';   #X-Coordinates
    y_vec=[0, 0, 0, 0, 0, 0];        #Y-Coordinates
end
