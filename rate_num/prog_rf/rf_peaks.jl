# This function finds the temporal values for starting and end of V values.
function compute_averge_V_trace(store_Varbar, fl)

    Tmin = netpar.TB * ceil(store_Varbar[1][1] / netpar.TB) - 100.0
    if Tmin < store_Varbar[1][1]
        println(fl.out, "Tmin=", Tmin, " < store_Varbar[1][1]=",
	  store_Varbar[1][1])
	exit(0)
    end
    Tmax = netpar.TB * floor(store_Varbar[end][1] / netpar.TB) - 100.0
    nTmin = ceil(Int64, (Tmin - runpar.epsilon) / runpar.deltat)
    nTmax = ceil(Int64, (Tmax - runpar.epsilon) / runpar.deltat)
    nTstart = ceil(Int64, (store_Varbar[1][1] - runpar.epsilon) / runpar.deltat)
    nTB = ceil(Int64, netpar.TB / runpar.deltat - runpar.epsilon)
    ncycles = div(nTmax - nTmin, nTB)
    printfmtln(fl.out, "Tmin={1:f} Tmax={2:f} nTmin={3:d} nTmax={4:d} " *
      "nTstart={5:d} nTB={6:d} ncycle={7:d}", Tmin, Tmax, nTmin, nTmax,
      nTstart, nTB, ncycles)

#   println(fl.out, "  ")
    VV_t_all = Array{Array{Float64, 1},1}(undef, nTB)
    [VV_t_all[iTB] = zeros(2) for iTB in 1:nTB]
    for iTB in 1:nTB
	iTindex = nTmin-nTstart+iTB
#	println("iTB=", iTB, " iTindex=", iTindex)
        VV_t_all[iTB][1] = store_Varbar[iTindex][1]
	VV_t_all[iTB][2] = 0.0
	for icycles in 1:ncycles
            VV_t_all[iTB][2] += store_Varbar[(icycles-1) * nTB + iTindex][2]
	end
	VV_t_all[iTB][2] /= ncycles
#	println(store_Varbar[iTindex][1], " ", store_Varbar[iTindex][2])
	printfmtln(fl.out, "{1:f} {2:f}", VV_t_all[iTB][1],
	  VV_t_all[iTB][2])
    end

    return VV_t_all, Tmin, Tmax
    
end #compute_averge_V_trace

# This function finds the peaks of V.
function find_V_peaks(VV_t_all, Tmin, Tmax, fl)

    length_half_chain = 4
    Vdiff = 2.0e-6
    println(fl.out, " ")

    val_max_ar = []
    val_min_ar = []
    
    for iTB in length_half_chain+1:length(VV_t_all)-length_half_chain
        Var = map(xx -> xx[2],
	  VV_t_all[iTB-length_half_chain:iTB+length_half_chain])
	  
        iTB_argmax = argmax(Var)
	Vmax = Var[iTB_argmax]
	
        iTB_argmin = argmin(Var)
	Vmin = Var[iTB_argmin]

	iTB_argmax == length_half_chain && VV_t_all[iTB][1] > Tmin + 100.0 &&
	  Vmax - Vmin > Vdiff ?
	  push!(val_max_ar, [VV_t_all[iTB][1], VV_t_all[iTB][2]]) : nothing
	  
	iTB_argmin == length_half_chain && Vmax - Vmin > Vdiff ?
	  push!(val_min_ar, [VV_t_all[iTB][1], VV_t_all[iTB][2]]) : nothing
    end

    val_ar_list = [val_max_ar, val_min_ar]
    for (i_ar_list, val_ar) in enumerate(val_ar_list)
        for rec in val_ar
            printfmtln(fl.out, "{1:f} {2:f} {3:d}", rec[1], rec[2], i_ar_list)
	end
	println(fl.out, "  ")
    end

    if length(val_min_ar) > 0 && length(val_max_ar) > 0
        if val_min_ar[1][1] > val_max_ar[1][1]
            println("val_min_ar[1][1]={1:f} > val_max_ar[1][1]={2:f}",
	      val_min_ar[1][1], val_max_ar[1][1])
            exit(0)
        end

        for iamp in 1:min(length(val_min_ar), length(val_max_ar))
            printfmtln(fl.wsk, "{1:d} {2:f} {3:f} {4:f}", iamp,
	      val_max_ar[iamp][2] - val_min_ar[iamp][2], val_max_ar[iamp][2],
	      val_min_ar[iamp][2])
        end
        flush(fl.wsk)
    
        for iamp in 1:length(val_max_ar)-1
            printfmtln(fl.tim, "{1:d} {2:f} {3:f} {4:f}", iamp,
	      val_max_ar[iamp+1][1] - val_max_ar[iamp][1],
	      val_max_ar[iamp+1][1], val_max_ar[iamp][1])
        end
        flush(fl.tim)
    end
    
end #find_V_peaks

# This function processed the values of V.
function processed_V_values(store_Varbar, fl)

    VV_t_all, Tmin, Tmax = compute_averge_V_trace(store_Varbar, fl)
    find_V_peaks(VV_t_all, Tmin, Tmax, fl)

end #processed_V_values
