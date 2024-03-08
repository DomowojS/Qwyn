#= 
Script to Excecute & Initialise 
=#
module Preprocessing_Input
export generateWF

function generateWF(prefix::AbstractString, folder::AbstractString)
    # Get a list of all files in the Inputfolder
    files = readdir(folder)
    # Filter out files that start with the given prefix
    script_files = filter(file -> occursin(r"^$prefix", file), files)
    # Execute each script
    for script_file in script_files
        include(joinpath(folder, script_file))
    end
end

generateWF("WF", "01_Input")

end