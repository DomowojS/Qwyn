#= 
Script to Excecute & Initialise 
=#
include("PKG_Manager.jl") #make neccecary packages available
module Preprocessing
export *

function generateWF(prefix::AbstractString, folder::AbstractString)
    # Get a list of all files in the Inputfolder
    files = readdir(folder)
    # Filter out files that start with the given prefix
    script_files = filter(file -> startswith(file, prefix), files)
    # Execute each script
    i::Int=0;
    #for script_file in script_files
        #i=i+1;
        #name="WF$i";
        #include(joinpath(folder, script_file))
        ## Create an instance of the struct
         #name= WF()
    #end
return i
end

end