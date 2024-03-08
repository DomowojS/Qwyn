#= Script for Package import
    This script imports all packages needed for the WF tool.
=#

# Working with directories
using Pkg
function PKG_Manager()
Pkg.add("FileIO")
end
using FileIO
