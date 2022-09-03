# This program finds the peaks of rf.col, folded within TB.

using Formatting
using Printf
using DataStructures

# This function reads the data in rf.col.* into the array store_Varbar
function read_col(store_Varbar, fl)

    for line in eachline(fl.col)
        line_split = parse.(Float64, split(line))
	push!(store_Varbar, [line_split[1], line_split[11]])
    end

#  for rec in store_Varbar
#      printfmtln(fl.out, "{1:f} {2:f}", rec[1], rec[2])
#  end
    
end #read_col

# main

dir_dat = pwd() * "/" #"../dat/"
include("rf_types.jl")
include("rf_peaks.jl")

if length(ARGS) >= 1
    suffix = ARGS[1]
else
    suffix = "a1"
end

fl = fl_st()
!isdir(suffix) ? mkdir(suffix) : nothing
fl.out = open(suffix * "/pk.out", "w")
fl.col = open(suffix * "/rf.col", "r")
fl.wsk = open(suffix * "/pk.wsk", "w")
fl.tim = open(suffix * "/pk.tim", "w")

include(dir_dat * suffix * "/rf_const.jl")

store_Varbar = Array{Array{Float64,1},1}(undef, 0) #2*netpar.Npop+1)

read_col(store_Varbar, fl)
processed_V_values(store_Varbar, fl)

close(fl.col)
close(fl.wsk)
close(fl.tim)
close(fl.out)
