# This file includes the functions used for analyzing bursting in vIRT neurons.

# This function finds the beginning and end times of population bursts.
function find_burst_initial_end_times(burst_init_end_times, ipop, rst_ar,
    time_between_bursts, runpar, fl)

    t_begin = runpar.Tall - runpar.tstat
    one_burst_time = [-1.0, -1.0]
    n_burst = 0

    for ifr in 1:length(rst_ar)
        if rst_ar[ifr].tspk - t_begin > time_between_bursts
	    if n_burst >= 1
	        one_burst_time[2] = t_begin
		to_put = zeros(2)
		[to_put[ii] = one_burst_time[ii] for ii in 1:2]
		push!(burst_init_end_times, to_put)
#		println(fl.out, "burst_init_end_times=", burst_init_end_times)
#		println(fl.out, n_burst,  " ", one_burst_time)
#		println("len=", length(burst_init_end_times))
	    end	    
            one_burst_time[1] = rst_ar[ifr].tspk
	    n_burst += 1
	    one_burst_time[1] = rst_ar[ifr].tspk
#	    println("n=", n_burst, " t=", t_begin, " ", rst_ar[ifr].tspk)
	end
	t_begin = rst_ar[ifr].tspk
    end
    println(fl.out, "ipop=", ipop, " ", n_burst,  " ", one_burst_time)
   

    if runpar.Tall - t_begin > time_between_bursts
    	if n_burst >= 1
	    one_burst_time[2] = t_begin
	    to_put = zeros(2)
            [to_put[ii] = one_burst_time[ii] for ii in 1:2]
            push!(burst_init_end_times, to_put)
	    println(fl.out, n_burst,  " ", one_burst_time)
	end	    
    end

    len_bur_time_ar = length(burst_init_end_times)
    printfmtln(fl.out, "length(burst_init_end_times)={1:d}",
      length(burst_init_end_times))
    
    if runpar.sm
        for irecord in 1:length(burst_init_end_times)
            record = burst_init_end_times[irecord]
            println(fl.out, irecord, " ", record)
	end
    end

    return len_bur_time_ar
end #find_burst_initial_end_times

# This function calculates the bursting time period.
function calculate_bursting_period(burst_init_end_times, runpar, fl)

    len_bur_ar = length(burst_init_end_times)
    if len_bur_ar >= 2
        Tper = (burst_init_end_times[end][2] - burst_init_end_times[1][2]) /
	  (len_bur_ar - 1)
    else
        Tper = -999.8
   end
   runpar.sm ? printfmtln(fl.out, "Tper={1:f}", Tper) : nothing

   return Tper
end #calculate bursting_period

# This function fits an envelope to the PSTH.
function generate_psth(rst_ar, Ncells, burstr, runpar, Npop_vIRT, fl)

    burstr.nbin =  floor(Int64, (runpar.tstat + runpar.epsilon) /
      runpar.dt_bin) 
    runpar.sm ? printfmtln(fl.out, "burstr.nbin={1:d}", burstr.nbin) : nothing
    psth_ar = Array{Array{Float64,1},1}(undef, Npop_vIRT)
    [psth_ar[ipop] = zeros(burstr.nbin+1) for ipop in 1:Npop_vIRT]

    for ipop in 1:Npop_vIRT
        for spk_tuple in rst_ar[ipop]
	    if spk_tuple.tspk > runpar.Tall - runpar.tstat - runpar.epsilon
	        tspk_for_bin = spk_tuple.tspk + runpar.tstat - runpar.Tall
	        itbin = floor(Int64, tspk_for_bin / runpar.dt_bin) + 1
	        deltat = tspk_for_bin - (itbin-1) * runpar.dt_bin
#	        println("t=", tspk_for_bin, " itbin=", itbin, " deltat=",
#                 deltat)
	        itbin > burstr.nbin ? (println("itbin=", itbin,
		  " > burstr.nbin=", burstr.nbin) ; exit(0)) : nothing
#	        println("itbin=", itbin)  
	        psth_ar[ipop][itbin]   += (runpar.dt_bin - deltat) /
		  runpar.dt_bin
	        psth_ar[ipop][itbin+1] +=           deltat  / runpar.dt_bin
	    end
	end
    end
    
#   [println("sum=", sum(psth_ar[ipop])) for ipop in 1:Npop_vIRT]
    for ipop in 1:Npop_vIRT
        for itbin in 1:burstr.nbin+1
	    psth_ar[ipop][itbin] *= 1.0e3 / (Ncells[ipop] * runpar.tstat)
	end
    end
#   [println("sum=", sum(psth_ar[ipop])) for ipop in 1:Npop_vIRT]

    prnt_bur = false
    if runpar.sm && prnt_bur
        for ipop in 1:Npop_vIRT
            for itbin in 1:burstr.nbin+1
	        printfmtln(fl.bur, "{1:f} {2:f} {3:d}",
		  (itbin-1.0) * runpar.dt_bin, psth_ar[ipop][itbin], itbin)
	    end
	    println(fl.bur, "  ")
        end
    end
   
    return psth_ar
end #generate_psth

# This function tests whether the point in the correlation array is a local
# extremum.
function is_it_exteremun(iar, cc_ar, n_test_extremum, imm, runpar)

    it_is_extremum = true
    index_list = collect(iar-n_test_extremum:iar-1)
    append!(index_list, iar+1:iar+n_test_extremum)
    for jar in index_list
#	    print("iar=", iar, " ", imm, " jar=", jar, " it=", it_is_extremum)
#	    println(" diff=", cc_ar[jar], " ", cc_ar[iar], " ",
#             cc_ar[jar] - cc_ar[iar], " xx=", imm * (cc_ar[jar] - cc_ar[iar]))
        if imm * (cc_ar[jar] - cc_ar[iar]) < 0.0
	    it_is_extremum = false
	end
    end

#    if iar >= 2 && cc_ar[iar] < runpar.epsilon &&
#      cc_ar[iar-1] < runpar.epsilon
#        it_is_extremum = false
#    end

#   println("it_is_extremum=", it_is_extremum)
    return it_is_extremum
end #is_it_exteremun

# This function checks if there are remaining indices of extrema.
function indices_remain(indices, indices_max)

    remain = false

    for ii in 1:length(indices)
        if indices[ii] < indices_max[ii]
	    remain = true
	end
    end

    return remain
end

# This function finds the index of the array with the smallest t (lag).
function find_index_of_smallest_t(indices, indices_max, ar_local_min_max)

    ii_t_min = 0
    for ii in 1:length(indices)
       if indices[ii] < indices_max[ii]
           ii_t_min = ii
	   break
       end
    end
    
    for ii in ii_t_min+1:length(indices)
        if (indices[ii] < indices_max[ii])
	    if ar_local_min_max[ii][indices[ii]+1] <
	       ar_local_min_max[ii_t_min][indices[ii_t_min]+1]
	         ii_t_min = ii
            end
	end
    end
    
    return ii_t_min
end #find_index_of_smallest_t

# This function combines the lists of local minima and maxima
function combine_ar_local(local_extremum_ar, ar_local_min_max, ipop, jpop,
    runpar, fl)

#   min, max, global_max
    indices_max = [length(ar_local_min_max[1]), length(ar_local_min_max[2]), 1]
    runpar.sm ? println(fl.out, "indices_max=", indices_max) : nothing
    indices = [0, 0, 0]

    while indices_remain(indices, indices_max)
        index_small_t = find_index_of_smallest_t(indices, indices_max,
	  ar_local_min_max)
	indices[index_small_t] += 1
	imm = index_small_t == 1 ? 1 : -1 #1 - min, -1 - max
	indices_here = indices[index_small_t]
	record = [imm, ar_local_min_max[index_small_t][indices_here]]
	push!(local_extremum_ar, record)
    end
    
    runpar.sm ? println(fl.out, "local_extremum_ar=", local_extremum_ar) :
      nothing

#   Consistency check
    for irecord in 2:length(local_extremum_ar)
        if local_extremum_ar[irecord-1][1] * local_extremum_ar[irecord][1] > 0
	    println("ipop=", ipop, " jpop=", jpop, " irecord=", irecord,
	      " local_extremum_ar=", local_extremum_ar[irecord-1],
	      local_extremum_ar[irecord])
	    exit(0)
	end
    end
    
end #combine_ar_local

# This function removes the spurious extrema that are "too local".
function remove_spurious_extrema(local_extremum_ar, cc_ar, i_argmax, ipop, jpop,
    runpar, fl)

    frac_height = 0.3

    mean_cc_ar = mean(cc_ar)
    std_cc_ar = Statistics.std(cc_ar)
    printfmtln(fl.out, "ipop={1:d} jpop={2:d} mean_cc_ar{3:f} std_cc_ar={4:f}",
      ipop, jpop, mean_cc_ar, std_cc_ar)

    find_spurious = false

    spurious_ar = []
    for irecord in 2:length(local_extremum_ar)
        cc1 = cc_ar[local_extremum_ar[irecord-1][2]]
	cc2 = cc_ar[local_extremum_ar[irecord][2]]
        if abs(cc2 - cc1) < frac_height * std_cc_ar
	    if find_spurious
	        push!(spurious_ar, irecord)
	    else
	        spurious_ar = [irecord-1, irecord]
            end
	    find_spurious = true
        else
	    if find_spurious
	        if isodd(length(spurious_ar))
                    it_in_spur = map(x -> local_extremum_ar[x][2], spurious_ar)
                    if i_argmax in it_in_spur
#		        println("i_argmax=", i_argmax, " it_in_spur=",
#			  it_in_spur)
		        for ii in spurious_ar
		            local_extremum_ar[ii][2] != i_argmax ? 
			      local_extremum_ar[ii][1] = 0 : nothing
			end
		    else
	                n_spur_half = div(length(spurious_ar), 2)
#		        local_extremum_ar[spurious_ar[n_spur_half+1]][1] *= -1
#		        loc_ext_n =
#			  local_extremum_ar[spurious_ar[n_spur_half+1]][2]
#		        loc_ext_m =
#			  local_extremum_ar[spurious_ar[n_spur_half]][2]
#		        loc_ext_p =
#			  local_extremum_ar[spurious_ar[n_spur_half+2]][2]
#		        print("itt=", local_extremum_ar[spurious_ar[n_spur_half+1]]
#		          , " ", local_extremum_ar[spurious_ar[n_spur_half]],
#		          " ", local_extremum_ar[spurious_ar[n_spur_half+2]],
#		          " ", loc_ext_n, " ", loc_ext_m, " ", loc_ext_p, " ",
#		          cc_ar[loc_ext_n], " ", cc_ar[loc_ext_m], " ",
#                         cc_ar[loc_ext_p])
#		        cc_ar[loc_ext_n] = 0.5 * (cc_ar[loc_ext_m] +
#		          cc_ar[loc_ext_p])
#                       println(" ", " cc=", cc_ar[loc_ext_n])
		        cc_in_spur = map(x -> cc_ar[local_extremum_ar[x][2]],
			  spurious_ar)
			
                        if mean(cc_in_spur) >= mean_cc_ar
		            i_cc_argmax = argmax(cc_in_spur)
			else
			    i_cc_argmax = argmin(cc_in_spur)
			end
                        list_to_remove = deleteat!(deepcopy(spurious_ar),
			  i_cc_argmax)
#		          n_spur_half+1)

                        for ii in list_to_remove
		            local_extremum_ar[ii][1] = 0
			end
		    end
		else
 		    println("i_argmax=", i_argmax, " spurious_ar=", spurious_ar)
		    if spurious_ar == [1, 2]
		        for ii in spurious_ar
		            local_extremum_ar[ii][1] = 0
			end
		    elseif spurious_ar == [length(local_extremum_ar)-1,
		      length(local_extremum_ar)]
		        for ii in spurious_ar
		            local_extremum_ar[ii][1] = 0
			end
		    elseif spurious_ar[1] > 1 &&
		      local_extremum_ar[spurious_ar[1]-1][2] == i_argmax
#		        println("loc_e=", local_extremum_ar[spurious_ar[1]-1],
#			  " i_argmax=", i_argmax)
                        for ii in spurious_ar
		            local_extremum_ar[ii][1] = 0
			end
	            elseif spurious_ar[end] < length(local_extremum_ar) &&
		      local_extremum_ar[spurious_ar[end]+1][2] == i_argmax
#		        println("loc_e=", local_extremum_ar[spurious_ar[end]+1],
#			  " i_argmax=", i_argmax)			  
                        for ii in spurious_ar
		            local_extremum_ar[ii][1] = 0
			end
	            elseif length(spurious_ar) == 2 &&
		      spurious_ar[2] - spurious_ar[1] == 1
		        for ii in spurious_ar
		            local_extremum_ar[ii][1] = 0
			end
                    elseif spurious_ar[end] - spurious_ar[1] ==
		      length(spurious_ar) - 1	
#	            elseif length(spurious_ar) == 4 &&
#		      spurious_ar[4] - spurious_ar[1] == 3
		        for ii in spurious_ar
		            local_extremum_ar[ii][1] = 0
			end
                    else
		        println("odd len(spurious_ar), no global maximum " *
			  "nearby. ipop=", ipop, " jpop=", jpop)
			exit(0)
	            end
		end
		spurious_ar = []
	        find_spurious = false
	    end
	end
    end
    
    if find_spurious
#	println("spar=", spurious_ar)
        for ii in spurious_ar
             local_extremum_ar[ii][1] = 0
	end
    end

    len_lo_ext = length(local_extremum_ar)
    for irecord in len_lo_ext:-1:1
        if local_extremum_ar[irecord][1] == 0
            deleteat!(local_extremum_ar, irecord)
	end
    end
    
end #remove_spurious_extrema

# This function finds the local maxima of the correlation function.
function find_local_maxima(locmax_ar, cc_ar, ipop, jpop, runpar, fl)

#   Find global maximum
    i_argmax = argmax(cc_ar)
    printfmtln(fl.out, "i_argmax={1:d} max={2:f}", i_argmax, cc_ar[i_argmax])
    
#   Find local maxima and minima
    n_test_extremum = 1
    ar_local_min = Array{Int64,1}(undef,0)
    ar_local_max = Array{Int64,1}(undef,0)
    ar_local_min_max = [ar_local_min, ar_local_max, [i_argmax]]

    for iminmax in 1:2
        list_indx = collect(n_test_extremum+1:i_argmax-1) # 19.10.2020 -2 -> -1
        append!(list_indx, i_argmax+1:2*runpar.lags-1)    # 19.10.2020 +2 -> +1

        for iar in list_indx
	    imm = 2*(iminmax%2)-1
            if is_it_exteremun(iar, cc_ar, n_test_extremum, imm, runpar)
	        push!(ar_local_min_max[iminmax], iar)
	    end
        end
        runpar.sm ? println(fl.out, "iminmax=", iminmax, " ar_local_min_max=",
	  ar_local_min_max[iminmax]) : nothing
    end

    local_extremum_ar = []
    combine_ar_local(local_extremum_ar, ar_local_min_max, ipop, jpop, runpar,
      fl)

    if runpar.sm 
        for record in local_extremum_ar
            time = (record[2] - 1 - runpar.lags) * runpar.dt_bin
	    cc = cc_ar[record[2]]
            printfmtln(fl.bur, "{1:f} {2:f} {3:d} {4:d} {5:d}", time,
	      1000.0 * cc, record[2], ipop, jpop)  
        end
        println(fl.bur, " ")
    end

    remove_spurious_extrema(local_extremum_ar, cc_ar, i_argmax, ipop, jpop,
      runpar, fl)

    if runpar.sm 
        for record in local_extremum_ar
            time = (record[2] - 1 - runpar.lags) * runpar.dt_bin
	    cc = cc_ar[record[2]]
            printfmtln(fl.bur, "{1:f} {2:f} {3:d} {4:d} {5:d}", time,
	      1000.0 * cc, record[2], ipop, jpop)  
        end
        println(fl.bur, " ")
    end

#   Find index of global maximum
    ind_global_max = 0
    for ii in 1:length(local_extremum_ar)
        if local_extremum_ar[ii][2] == i_argmax
	    ind_global_max = ii
	    break
	end
    end
    
    if ind_global_max == 0
        println("ipop=", ipop, " jpop=", jpop, " ind_global_max=",
	  ind_global_max, " i_argmax=", i_argmax)
#	exit(0)
    end
#   println("ind_global_max=", ind_global_max, " it=",
#     local_extremum_ar[ind_global_max][2])

    return ind_global_max, local_extremum_ar
end #function find_local_maxima

# This function finds Tper from the cross-correlation.
function find_Tper_from_cc(ind_global_max, local_extremum_ar, cc_ar, ipop, jpop,
    runpar, fl)

    lag_used_t_min = max(ind_global_max - runpar.lag_for_Tper, 1)
    lag_used_t_max = min(ind_global_max + runpar.lag_for_Tper,
      length(local_extremum_ar))
    if lag_used_t_max > lag_used_t_min
        t_small = local_extremum_ar[lag_used_t_min][2] * runpar.dt_bin
        t_large = local_extremum_ar[lag_used_t_max][2] * runpar.dt_bin
        Tper = 2.0 * (t_large - t_small) / (lag_used_t_max - lag_used_t_min)
	if runpar.sm
            print(fl.out, " ind_global_max=", ind_global_max, " lag_used_t_min=",
              lag_used_t_min, " lag_used_t_max=", lag_used_t_max)
            println(fl.out, " t_small=", t_small, " t_large=", t_large,
	      " Tper=", Tper)
	end
    else
        Tper = -999.6
   end
   
   return Tper
end #find_Tper_from_cc

# This function computes the cross-correlation of the psth.
function burst_cross_cor(psth_ar, Ncells, burstr, runpar, Npop_vIRT, fl)

    cc_mat_ar = Array{Array{Array{Float64,1},1},1}(undef, Npop_vIRT)
    locmax_mat_ar = Array{Array{Array{Int64,1},1},1}(undef, Npop_vIRT)
    for ipop in 1: Npop_vIRT
        cc_mat_ar[ipop] = Array{Array{Float64,1},1}(undef, Npop_vIRT)
	locmax_mat_ar[ipop] = Array{Array{Int64,1},1}(undef, Npop_vIRT)
        for jpop in 1:ipop
            cc_mat_ar[ipop][jpop] = Array{Float64,1}(undef, 2*runpar.lags+1)
            locmax_mat_ar[ipop][jpop] = Array{Int64,1}(undef, 2*runpar.lags+1)
        end
    end

    tvec = zeros(2*runpar.lags+1)
#   Computing cross-corralations
    for ipop in 1:Npop_vIRT
        for jpop in 1:ipop
            cc_mat_ar[ipop][jpop] = crosscov(psth_ar[ipop], psth_ar[jpop],
	      -runpar.lags:runpar.lags, demean=false)
	      
	    for itbin in  1:2*runpar.lags+1
                time = (itbin - 1 - runpar.lags) * runpar.dt_bin
		tvec[itbin] = time
		
	        cc_mat_ar[ipop][jpop][itbin] *= runpar.tstat /
		  (runpar.tstat - abs(itbin - runpar.lags -  1))
	    end

            if runpar.sm 
                for itbin in 1:2*runpar.lags+1
                    time = (itbin - 1 - runpar.lags) * runpar.dt_bin
                    printfmtln(fl.bur, "{1:f} {2:f} {3:d} {4:d} {5:d}", time,
		      1.0e6 * cc_mat_ar[ipop][jpop][itbin], itbin, ipop, jpop)
	        end
                printfmtln(fl.bur, "  ")
	    end

#           Smoothing cross-correlations
            lambda_spline = 1.0e2
            spl = SmoothingSplines.fit(SmoothingSpline, tvec,
	      cc_mat_ar[ipop][jpop], lambda_spline)
            xsmo = SmoothingSplines.predict(spl) # fitted vector

            if runpar.sm 
                for itbin in 1:2*runpar.lags+1
                    time = (itbin - 1 - runpar.lags) * runpar.dt_bin
                    printfmtln(fl.bur, "{1:f} {2:f} {3:d} {4:d} {5:d}", time,
		      1.0e6 * xsmo[itbin], itbin, ipop, jpop)
	        end
                printfmtln(fl.bur, "  ")
	    end

            cc_mat_ar[ipop][jpop] = deepcopy(xsmo)

	end
    end

#   Find local maxima
    n_Tper_from_cc_av = 0
    Tper_from_cc_av = 0.0
    
    for ipop in 1:Npop_vIRT
        for jpop in 1:ipop
#if ipop == 2 && jpop == 1
	    ind_global_max, local_extremum_ar = find_local_maxima(
	      locmax_mat_ar[ipop][jpop], cc_mat_ar[ipop][jpop], ipop, jpop,
	      runpar, fl)
	    if ind_global_max != 0
                Tper_from_cc = find_Tper_from_cc(ind_global_max,
		  local_extremum_ar, cc_mat_ar[ipop][jpop], ipop, jpop, runpar,
		  fl)
	        printfmtln("ipop={1:d} jpop={2:d} Tper_from_cc={3:f}", ipop,
		  jpop, Tper_from_cc)
	        n_Tper_from_cc_av += 1
	        Tper_from_cc_av += Tper_from_cc
	    end
#end	    
	end
    end
    
    if n_Tper_from_cc_av > 0
        Tper_from_cc_av /= n_Tper_from_cc_av
    else
        Tper_from_cc_av = -999.5
    end
    runpar.sm ? println(fl.out, "Tper_from_cc_av=", Tper_from_cc_av) : nothing

    return Tper_from_cc_av
end #burst_cross_cor

# This function computes CV and CV2 for tonic firing neuronal populations.
function compute_CV_tonic(rst_ar, ipop, Ncells, mode, Tper_from_cc_av, runpar,
    fl)

    frac_tper_for_CV = 0.4
    nspk_ion = zeros(Int64, Ncells)
    n_deltat_spk_ion_CV = zeros(Int64, Ncells)
    n_deltat_spk_pair_CVt = zeros(Int64, Ncells)
    tspk_old = zeros(Ncells)
    delta_tspk_old = zeros(Ncells)
    isi_av = zeros(Ncells)
    isi_avt = zeros(Ncells)
    CV_ion = zeros(Ncells)
    CVt_ion = zeros(Ncells)

#   Skipping over the spikes before the tstat time interval.
    t_begin = runpar.Tall - runpar.tstat
    nspk_beg = 1
    while nspk_beg < length(rst_ar)
        rst_ar[nspk_beg].tspk < t_begin ? nspk_beg += 1 : break
    end
    
#   println("nspk_beg=", nspk_beg, " t_begin=", t_begin, " tspk=",
#     rst_ar[nspk_beg].tspk)
    
    for ispk in nspk_beg:length(rst_ar)
        ion_pop = rst_ar[ispk].ion_pop
        tspk = rst_ar[ispk].tspk
        nspk_ion[ion_pop] += 1

#       CV
	if nspk_ion[ion_pop] >= 2
	    delta_tspk = tspk - tspk_old[ion_pop]
#	    ipop == 1 && ion_pop == 1 ? println("\ndelta_tspk=", delta_tspk) : nothing
	    
	    if mode != 2 || delta_tspk < frac_tper_for_CV * Tper_from_cc_av
#	        println(" mode=", mode, " ion_pop=", " t=", tspk_old[ion_pop],
#                 tspk, " delta_tspk=", delta_tspk, " fr=", frac_tper_for_CV,
#                 " T=", Tper_from_cc_av, " fr*T=", frac_tper_for_CV *
#                 Tper_from_cc_av)
	        n_deltat_spk_ion_CV[ion_pop] += 1
	        isi_av[ion_pop] += delta_tspk
                isi_avt[ion_pop] += delta_tspk * delta_tspk
#		ipop == 1 && ion_pop == 1 ? println("C ",
#		  n_deltat_spk_ion_CV[ion_pop], " ", tspk_old[ion_pop], " ",
#		  tspk, " ", delta_tspk, " ", isi_av[ion_pop], " ",
#		  isi_avt[ion_pop]) : nothing

#               CVt
	        if (mode != 2 && nspk_ion[ion_pop] >= 3) ||
	          (n_deltat_spk_ion_CV[ion_pop] >= 2 &&
		  delta_tspk_old[ion_pop] < frac_tper_for_CV * Tper_from_cc_av)
		    n_deltat_spk_pair_CVt[ion_pop] += 1
                    CVt_now = 2.0* abs(delta_tspk - delta_tspk_old[ion_pop]) /
	              (delta_tspk + delta_tspk_old[ion_pop])
	            CVt_ion[ion_pop] += CVt_now
		    
#		    ipop == 1 && ion_pop == 1 ? println("t ",
#		      n_deltat_spk_pair_CVt[ion_pop], " dtsp=",
#		      delta_tspk_old[ion_pop], " ", delta_tspk, " CVt_now=",
#		      CVt_now, " CVt_ion=", CVt_ion[ion_pop]) : nothing	    
	        end
	    else
#	        println(" mode=", mode, "ion_pop=", " t=", tspk_old[ion_pop], tspk, " delta_tspk=", delta_tspk)
	    end
	    
            delta_tspk_old[ion_pop] = delta_tspk
#	    ipop == 1 && ion_pop == 1 ? println("dtsp=",
#	      delta_tspk_old[ion_pop]) : nothing
	end
	
        tspk_old[ion_pop] = tspk
#	ipop == 1 && ion_pop == 1 ? println("tsp=", tspk_old[ion_pop]) : nothing
    end

#   Computing statistics of CV and CVt.
    for ion_pop in 1:Ncells
        if n_deltat_spk_ion_CV[ion_pop] <= 1 #nspk_ion[ion_pop] <= 2
	    CV_ion[ion_pop] = -999.9
	    CVt_ion[ion_pop] = -999.8
	else
	    isi_av[ion_pop] /=  n_deltat_spk_ion_CV[ion_pop] #nspk_ion[ion_pop] - 1
	    isi_avt[ion_pop] /=  n_deltat_spk_ion_CV[ion_pop] #nspk_ion[ion_pop]- 1
	    diff = isi_avt[ion_pop] - isi_av[ion_pop] * isi_av[ion_pop]
	    if diff > 0
	        sd = sqrt(diff)
	    elseif diff > -runpar.epsilon
	        sd = 0.0
	    else
	        println("diff=", diff, " <= 0!")
		exit(0)
	    end
    	    CV_ion[ion_pop] = sd / isi_av[ion_pop]
	    
      	    if n_deltat_spk_pair_CVt[ion_pop] <= 1
	    # n_deltat_spk_ion_CV[ion_pop] <= 2 #nspk_ion[ion_pop] <= 3
	        CVt_ion[ion_pop] = -999.85
	    else
                CVt_ion[ion_pop] /= n_deltat_spk_pair_CVt[ion_pop]
#		n_deltat_spk_ion_CV[ion_pop] - 1 # nspk_ion[ion_pop] - 2
	    end
        end
    end

    if runpar.sm
        println(fl.out, "  ion    CV      CVt   isiav  n_deltat_spk_ion_CV  n_deltat_spk_pair_CVt")
        for ion_pop in 1:Ncells
            printfmtln(fl.out, "{1:d} {2:f} {3:f} {4:f} {5:d} {6:d}", ion_pop,
	      CV_ion[ion_pop], CVt_ion[ion_pop], isi_av[ion_pop],
	      n_deltat_spk_ion_CV[ion_pop], n_deltat_spk_pair_CVt[ion_pop])
	end
    end
    
#   Compute population-averaged values
    n_silent = 0
    n_no_CVt = 0
    CV_pop_av  = 0.0
    CVt_pop_av = 0.0
	
    for ion_pop in 1:Ncells
        if n_deltat_spk_ion_CV[ion_pop] <= 2 #nspk_ion[ion_pop] <= 3
            n_silent += 1
        else
            CV_pop_av += CV_ion[ion_pop]
	    if n_deltat_spk_pair_CVt[ion_pop] <= 2
	        n_no_CVt += 1
	    else
                CVt_pop_av += CVt_ion[ion_pop]
	    end
        end
    end

    if n_silent == Ncells
#       all silent
        CV_pop_av  = -998.0
    else
        CV_pop_av  /= Ncells - n_silent
    end

    if n_no_CVt == Ncells
#       no CVt defined
        CVt_pop_av = -998.1
    else
        CVt_pop_av /= Ncells - n_no_CVt
    end

#	av[ipop].n_silent   = ISIav.n_silent[ipop]
#	av[ipop].CV_pop_av  = ISIav.CV_pop_av[ipop]
#	av[ipop].CVt_pop_av = ISIav.CVt_pop_av[ipop]
#	av[ipop].T_pop_av   = ISIav.T_pop_av[ipop]
    
    printfmtln(fl.out, "n_silent={1:d} CV_pop_av={2:f} CV_pop_avt={3:f} ",
      n_silent, CV_pop_av, CVt_pop_av)

    return CV_pop_av, CVt_pop_av #, n_silent
end #compute_CV_tonic

# This function finds the maximal ISI for each neuron seperately, and then
# the population-statistics of these ISIs.
function find_max_isi_for_each_neuron(isi_max_neuron, rst_ar, Ncells, ipop, fl)

    t_spk_ion_ar = -ones(Ncells)

    nrst = 1
    while nrst < length(rst_ar) &&
      rst_ar[nrst].tspk < runpar.Tall - runpar.tstat - runpar.epsilon
        nrst += nrst
    end
#   println("nrst=", nrst, " tspk=", rst_ar[nrst].tspk)
    
    for irst in nrst:length(rst_ar)
        ion_pop = rst_ar[irst].ion_pop
	tspk =  rst_ar[irst].tspk
	
	if t_spk_ion_ar[ion_pop] > 0.0
	    delta_t_spk = tspk - t_spk_ion_ar[ion_pop]
	    delta_t_spk > isi_max_neuron[ion_pop] ? isi_max_neuron[ion_pop] =
	      delta_t_spk : nothing
	end
	
	t_spk_ion_ar[ion_pop] = tspk
    end

    positive_isi_max_neuron = filter(x -> x > 0.0, isi_max_neuron)
    if (length(positive_isi_max_neuron)) >= 1
        isi_max_neuron_median = median(positive_isi_max_neuron)
    else
        isi_max_neuron_median = -995.4
    end
    
    if runpar.sm
        printfmtln(fl.out, "isi_max_neuron ipop={1:d}  median={2:f}", ipop,
	  isi_max_neuron_median)
	for ion_pop in 1:Ncells
	    printfmtln(fl.out, "{1:d} {2:f}", ion_pop, isi_max_neuron[ion_pop])
	end
    end

    return isi_max_neuron_median
end #find_max_isi_for_each_neuron

# This function carries out the analysis of bursting.
function burst_analyze(burstr, rst_ar, Ncells, runpar, av, fl)

    Npop = length(Ncells)
    Npop_vIRT = 2
    ratio_max_av = 0.5
    ratio_max_av_mode = 2.0

    isi_max = zeros(Npop_vIRT)
    fr_av = zeros(Npop)
    isi_av = zeros(Npop_vIRT)
    time_between_bursts = zeros(Npop_vIRT)
    burst_init_end_times = Array{Array{Array{Float64,1},1},1}(undef,Npop_vIRT)
    [burst_init_end_times[ipop] = Array{Array{Float64,1}}(undef,0) for
      ipop in 1:Npop_vIRT]
    Tper_ar = zeros(Npop_vIRT)
    isi_max_neuron = Array{Array{Float64,1}}(undef,Npop_vIRT)
    [isi_max_neuron[ipop] = zeros(Ncells[ipop]) for ipop in 1:Npop_vIRT]

    for ipop in 1:Npop #_vIRT 
        ispk = 1
        while ispk <= length(rst_ar[ipop])
	    if rst_ar[ipop][ispk].tspk > runpar.Tall - runpar.tstat -
	      runpar.epsilon
	        break
	    end
	    ispk += 1
	end
	
        fr_av[ipop] = 1.0 * (length(rst_ar[ipop]) - ispk + 1) / runpar.tstat /
	  Ncells[ipop]

        if ipop <= Npop_vIRT 
	    if fr_av[ipop] > runpar.epsilon
                isi_av[ipop] = fr_av[ipop] > 0 ?  1.0 / fr_av[ipop] : -999.9
	        len_ar = length(rst_ar[ipop])
	        diff_ar = zeros(len_ar-1)
	        [diff_ar[it] = rst_ar[ipop][it+1].tspk - rst_ar[ipop][it].tspk
		  for it in 1:len_ar-1]
	      
#	        isi_max[ipop] = extrema(diff_ar)[2]
                if fr_av[ipop] < runpar.epsilon
	            isi_max[ipop] = -995.5
	        else
	            isi_max[ipop] = find_max_isi_for_each_neuron(
		      isi_max_neuron[ipop], rst_ar[ipop], Ncells[ipop], ipop,
		      fl)
		end
	    end
	end
    end

    printfmtln(fl.out, "fr_av={1:f} {2:f} {3:f}", fr_av[1], fr_av[2], fr_av[3])
    printfmtln(fl.out, "isi_av={1:f} {2:f}", isi_av[1], isi_av[2])
    printfmtln(fl.out, "isi_max={1:f} {2:f}", isi_max[1], isi_max[2])

    burstr.Mfr = zeros(Npop)
    [burstr.Mfr[ipop] = fr_av[ipop] for ipop in 1:Npop]
    
    if fr_av[1] < runpar.epsilon && fr_av[2] < runpar.epsilon
#       Silent network.
        burstr.mode = 0
        burstr.Tper = -999.4
        burstr.Tper_from_cc_av = -999.4
    elseif (fr_av[1] < runpar.epsilon && fr_av[2] > runpar.epsilon) ||
           (fr_av[1] > runpar.epsilon && fr_av[2] < runpar.epsilon)
#       One population is silent, one is active.	   
        burstr.mode = 3
        burstr.Tper = -999.5
        burstr.Tper_from_cc_av = -999.5
    elseif isi_max[1] < ratio_max_av_mode * isi_av[1] ||
           isi_max[2] < ratio_max_av_mode * isi_av[2]
#       The two populations fire, but do not burst.
        burstr.mode = 1
        burstr.Tper = -999.7
        burstr.Tper_from_cc_av = -999.7
    else
#       The two populations burst.
        burstr.mode = 2

        time_between_bursts = ratio_max_av .* isi_max +
          (1.0 - ratio_max_av) * isi_av
      
        printfmtln(fl.out, "time_between_bursts={1:f} {2:f}",
          time_between_bursts[1], time_between_bursts[2])

        for ipop in 1:Npop_vIRT
            len_bur_time_ar = find_burst_initial_end_times(
	      burst_init_end_times[ipop], ipop,
	      rst_ar[ipop], time_between_bursts[ipop], runpar, fl)
	      
#	    if len_bur_time_ar <= 2
#	         burstr.mode += 10
#	    end
	    
            Tper_ar[ipop] = calculate_bursting_period(
	      burst_init_end_times[ipop], runpar, fl)
        end

        psth_ar = generate_psth(rst_ar, Ncells, burstr, runpar, Npop_vIRT,
	  fl)
#       fit_burst_profile(psth_ar, Ncells, burstr, runpar, Npop_vIRT, fl)

        Tper_from_cc_av = burst_cross_cor(psth_ar, Ncells, burstr, runpar,
	  Npop_vIRT, fl)
	
        burstr.Tper = (Tper_ar[1] + Tper_ar[2]) / 2.0
	burstr.Tper_from_cc_av = Tper_from_cc_av
    end
 
    printfmt(fl.out, "mode={1:d} Tper={2:f} Mfr=", burstr.mode, burstr.Tper)
    [printfmt(fl.out, "{1:f} ", burstr.Mfr[ipop]) for ipop in 1:Npop_vIRT]
    printfmt(fl.out, "\n")

    av[1].mode = burstr.mode
    av[1].Tper = burstr.Tper
    [av[ipop].Mfr = burstr.Mfr[ipop] for ipop in 1:Npop]
    av[1].Tper_from_cc_av = burstr.Tper_from_cc_av
	println(fl.out, "burstr.Tper_from_cc_av=", burstr.Tper_from_cc_av,
	" av=", av[1].Tper_from_cc_av)

    if burstr.mode == 0
        for ipop in 1:Npop_vIRT
            av[ipop].CV_bur_pop_av = -996.1
	    av[ipop].CVt_bur_pop_av = -996.1
	end
    elseif burstr.mode == 1
        for ipop in 1:Npop_vIRT
	    av[ipop].CV_bur_pop_av, av[ipop].CVt_bur_pop_av = compute_CV_tonic(
	      rst_ar[ipop], ipop, Ncells[ipop], burstr.mode,
	      burstr.Tper_from_cc_av, runpar, fl)
	end
    elseif burstr.mode == 2
        for ipop in 1:Npop_vIRT
	    av[ipop].CV_bur_pop_av, av[ipop].CVt_bur_pop_av = compute_CV_tonic(
	      rst_ar[ipop], ipop, Ncells[ipop], burstr.mode,
	      burstr.Tper_from_cc_av, runpar, fl)
	end
    elseif burstr.mode == 3
        if fr_av[1] < runpar.epsilon
	    av[1].CV_bur_pop_av = -996.2
	    av[1].CVt_bur_pop_av = -996.2
	    av[2].CV_bur_pop_av, av[2].CVt_bur_pop_av = compute_CV_tonic(
	      rst_ar[2], 2, Ncells[2], burstr.mode, burstr.Tper_from_cc_av,
	      runpar, fl)
	else
	    av[2].CV_bur_pop_av = -996.2
	    av[2].CVt_bur_pop_av = -996.2
	    av[1].CV_bur_pop_av, av[1].CVt_bur_pop_av = compute_CV_tonic(
	      rst_ar[1], 1, Ncells[1], burstr.mode, burstr.Tper_from_cc_av,
	      runpar, fl)
	end
    end
    
end # burst_analyze
