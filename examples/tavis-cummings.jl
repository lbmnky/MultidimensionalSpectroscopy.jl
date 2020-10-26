#!/usr/bin/julia
# using PyPlot
# using QuantumOptics

# make sure to set script directory as pwd()
cd(@__DIR__)

# include my custom cmds module
if Sys.iswindows()
    include("..\\cmds.jl")
    fn = "01_Output"
else
    include("../cmds.jl")
    fn = "01_Output"
end

pygui(true)

calc_2d = true

cmp = cmds.create_colormap("bright");
