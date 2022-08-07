# This file includes the functions used for computing the hights of peaks
# of whisker protractions after bursting of a preBot neuron.

using SmoothingSplines

# Sigmoid function
sigmoid(xx, nn) = xx^nn / (1 + xx^nn)

# This function computes the whisker angle using a linear model.
function compute_whisker_angle_linear(theta_t, peakpar, rst_ar, partohop, fl)
    tau1 = 5.0
    tau2 = 30.0
    
    time = 0.0
    ex1 = 0.0
    ex2 = 0.0
    ex1in = 0.0
    ex2in = 0.0

#   for line in eachline(fl.rst)
#       line_split = split(line, " ")
#       tspk = parse(Float64, line_split[1])
#       ipop = parse(Int64, line_split[3])

#   println("DD length_rst_ar=", length(rst_ar))
    for spike_tuple in rst_ar[partohop.ipop_cal_wsk]
        tspk = spike_tuple.tspk
#       peakpar.sm ? println(fl.out, tspk, " ", ipop,  " ", time) : nothing
    
        if tspk > time

            if tspk - time > peakpar.deltat
		    ntime = floor(Int64, (tspk - time) / peakpar.deltat)
		    if peakpar.sm
		        printfmtln(fl.out, "time={1:f} tspk={2:f} ntime={3:d}",
			  tspk, time, ntime)
	            end
		    for itime in 1:ntime
		        ex1in = ex1 * exp(-itime * peakpar.deltat / tau1)
			ex2in = ex2 * exp(-itime * peakpar.deltat / tau2)
			exall = (ex2in - ex1in) / (tau2 - tau1)
			tnow = time + itime * peakpar.deltat
			push!(theta_t, [tnow, exall])
		    end
		end

            ex1 *= exp(-(tspk - time) / tau1) 
            ex2 *= exp(-(tspk - time) / tau2) 
            ex1 += 1
    	    ex2 += 1
        elseif tspk - time > -peakpar.deltat_spk
            ex1 += 1
            ex2 += 1
        else
            println("tspk-time=", tspk-time)
        	exit(0)
        end
    
        exall = (ex2 - ex1) / (tau2 - tau1)
	push!(theta_t, [tspk, exall])
        time = tspk
    end

#   time from the last spike to the end of the time interval
    ntime = floor(Int64, (partohop.Tall - time) / peakpar.deltat)
#   println("DD2 Tall=", partohop.Tall, " time=", time, " ntime=", ntime)
    for itime in 1:ntime
        ex1in = ex1 * exp(-itime * peakpar.deltat / tau1)
    	ex2in = ex2 * exp(-itime * peakpar.deltat / tau2)
    	exall = (ex2in - ex1in) / (tau2 - tau1)
    	tnow = time + itime * peakpar.deltat
    	push!(theta_t, [tnow, exall])
    end
    
    ex1 *= exp(-(partohop.Tall - time) / tau1)
    ex2 *= exp(-(partohop.Tall - time) / tau2) 
    exall = (ex2 - ex1) / (tau2 - tau1)
    push!(theta_t, [partohop.Tall, exall])

end #compute_whisker_angle_linear

# This function writes the whisking trajectory theta(t) on fl.wsk .
function write_wsk(theta_t, peakpar, partohop, fl)

    if peakpar.sm && runpar.open_file_for_par
        for theta_t_record in theta_t
            printfmtln(fl.wsk, "{1:f} {2:f}", theta_t_record[1],
	      theta_t_record[2])
	end
    end
end #write_wsk 

# This function computes the average and sd of theta(t).
function compute_statistics_wsk(theta_t, peakpar, partohop, awsval, fl)

    t_cal = 0.0
    sum = 0.0
    sum_2 = 0.0
    itheta = 1

    while itheta <=length(theta_t)
       if theta_t[itheta][1] > partohop.Tall - partohop.tstat
           if itheta > 1
	       del_t_cal = theta_t[itheta][1] - partohop.Tall + partohop.tstat
	       val_former = lininter(theta_t[itheta-1][1], theta_t[itheta][1],
	         partohop.Tall - partohop.tstat,
		 theta_t[itheta-1][2], theta_t[itheta][2])
	       val = theta_t[itheta][2] + val_former
	       t_cal += del_t_cal
	       sum = val * del_t_cal / 2.0
	       sum_2 = val * val * del_t_cal / 2.0
	       itheta += 1
	       break
	   end
       end
       itheta += 1
    end

    while itheta < length(theta_t)
        itheta += 1
	del_t_cal = theta_t[itheta][1] - theta_t[itheta-1][1]
	val= (theta_t[itheta][2] + theta_t[itheta-1][2]) / 2.0
	t_cal += del_t_cal
	sum += val * del_t_cal
	sum_2 += val * val * del_t_cal
    end
    
    if t_cal > 0.0
        whisk_set = sum / t_cal
	diff = sum_2 / t_cal - whisk_set * whisk_set
	if diff > 0.0
	    whisk_amp = sqrt(diff)
	elseif diff > -partohop.epsilon
	    whisk_amp = 0.0
	else
	    whisk_amp = -998.7
	end
    else
        whisk_set = 0.0
	whisk_amp = 0.0
    end

    printfmtln(fl.out, "whisk_set={1:f} whisk_amp={2:f}", whisk_set, whisk_amp)
    awsval.whisk_set = whisk_set
    awsval.whisk_amp = whisk_amp
#   printfmt(fl.aws, "{1:f} {2:f}", whisk_set, whisk_amp)

    return whisk_set, whisk_amp
end #compute_statistics_wsk


# This function find the new index in the array of spikes of FN neurons
# for which the spike_time <= time
function find_new_index(rst_ar, index_FN_spike, time,
	  partohop, fl)

    FN_spikes_array = [ ]
    new_index_FN_spike = index_FN_spike

#   short-circut evaluation of the condition in the while loop.
    while new_index_FN_spike < length(rst_ar[partohop.ipop_cal_wsk]) &&
      rst_ar[partohop.ipop_cal_wsk][new_index_FN_spike+1].tspk <= time
        new_index_FN_spike += 1
	
	spike_tuple = rst_ar[partohop.ipop_cal_wsk][new_index_FN_spike]
	push!(FN_spikes_array, [spike_tuple.tspk, spike_tuple.ion_pop])
    end

    return (new_index_FN_spike, FN_spikes_array)
end


# This function computes the whisker angle using a nonlinear model for which
# the force is a fraction of powers of [Ca] (Simony et al. 2010).
function compute_whisker_angle_power(theta_t, peakpar, rst_ar, partohop, fl)

#   r0 = 1.9 * 0.1
#   taur = 5.0
#   tauc = 20.0
#   tauw = 20.0
#   AA = 100.0

    r0   = partohop.whisker.r0
    taur = partohop.whisker.taur
    tauc = partohop.whisker.tauc
    tauw = partohop.whisker.tauw
    AA   = partohop.whisker.AA

    rtrc = r0 * tauc / (tauc - taur)
    
    time = 0.0

    tspk_last = zeros(1, partohop.nonF)
    Ck = zeros(1, partohop.nonF)
    Ca = zeros(1, partohop.nonF)
    force_amp = zeros(1, partohop.nonF)
    index_FN_spike = 0
    theta_whisker = 0.0

    push!(theta_t, [0.0, 0.0])

    [Ck[ion] = 0.0 for ion in 1:partohop.nonF]

    while time < partohop.Tall - partohop.epsilon
        time += peakpar.deltat

#       Finding FN spikes in this time step.
        (new_index_FN_spike, FN_spikes_array) = find_new_index(rst_ar,
	  index_FN_spike, time, partohop, fl)
#       println(fl.out, time, " ", index_FN_spike, " ", new_index_FN_spike,
#       " ", index_FN_spike, " ", length(rst_ar))
#       println(fl.out, "FN_spikes_array=", FN_spikes_array)
	index_FN_spike = new_index_FN_spike

        ion_fire_set = Set(map(x -> Int64(round(x[2])), FN_spikes_array))
#	printfmt("time={1:f} l={2:d}", time, length(ion_fire_set))
#	println(" ion_fire_set=", ion_fire_set)

##      printfmt(fl.out, "{1:f}", time)
        for ion in 1:partohop.nonF
	    if ion in ion_fire_set
		tspk_new = filter(x -> Int64(round(x[2])) == ion,
		  FN_spikes_array)[1][1]
#		println("time=", time, " ion=", ion, " ts=", tspk_new)
		
	        if tspk_last[ion] > partohop.epsilon
		    exr = exp(-(tspk_new - tspk_last[ion])/taur)
		    exc = exp(-(tspk_new - tspk_last[ion])/tauc)
		    Ck[ion] = Ck[ion] * exc + rtrc * (exc - exr)
		end
		tspk_last[ion] = tspk_new 
		exr = exp(-(time - tspk_last[ion])/taur)
		exc = exp(-(time - tspk_last[ion])/tauc)
		Ca[ion] = rtrc * (exc - exr) + Ck[ion] * exc
	    else
	        if tspk_last[ion] > partohop.epsilon
		    exr = exp(-(time - tspk_last[ion])/taur)
		    exc = exp(-(time - tspk_last[ion])/tauc)
		    Ca[ion] = rtrc * (exc - exr) + Ck[ion] * exc
		end
	    end
	    force_amp[ion] = AA * sigmoid(Ca[ion], 4.0)
#	    printfmt(fl.out, " {1:f}", force_amp[ion]) #Ca[ion])
	end
	force_amp_all = sum(force_amp) / partohop.nonF
	theta_whisker += peakpar.deltat * (-theta_whisker / tauw +
	  force_amp_all)
#	printfmt(fl.out, "{1:f} {2:f} {3:f} {4:f} {5:f}\n", time,
#	  theta_whisker, force_amp_all, force_amp[1], Ca[1])
	push!(theta_t, [time, theta_whisker])
    end

#   println(fl.out, "  ")
    
end #compute_whisker_angle_power

# This function finds the starting times of preBot synaptic strength, assuming
# external preBot input.
#function find_time_external_prebot(prebot_burst_start, peakpar, partohop, fl)

#    for iper in 0:div(partohop.Tall, partohop.Tper_syn_ext)
#         T_start_per = iper * partohop.Tper_syn_ext
#         push!(prebot_burst_start, T_start_per)
#    end

#    peakpar.sm ? println(fl.out, "prebot_burst_start=", prebot_burst_start) :
#      nothing
    
#end # find_time_external_prebot

# This function computes the time course of the whisking angle, writes it
# on irt.wsk or hop.wsk, and compute the set point and the average whisking
# amplitude.
function compute_write_whisking_angle(theta_t, peakpar, rst_ar, partohop,
  awsval, fl)

    if partohop.angle_cal == 'l'
        compute_whisker_angle_linear(theta_t, peakpar, rst_ar, partohop, fl)
    elseif partohop.angle_cal == 'p'
        compute_whisker_angle_power(theta_t, peakpar, rst_ar, partohop, fl)
    else
        println("angle_cal=", partohop.angle_cal, " should be l or p!")
	exit(0)
    end

    write_wsk(theta_t, peakpar, partohop, fl)
    whisk_set, whisk_amp = compute_statistics_wsk(theta_t, peakpar, partohop,
      awsval, fl)

    return whisk_set, whisk_amp
end #compute_write_whisking_angle

#----------------------------------------------------

# This function takes the whisking trace in theta_t and extracts tvec (time) and
# xvec (theta).
function extract_smooth_whisking_trace(theta_t, peakpar, partohop, fl)

    lambda_spline = 1.0e4 # 0.2e-6

    tstart = partohop.Tall - partohop.tstat - partohop.epsilon
    tvec = map(x -> x[1], filter(x -> x[1] > tstart, theta_t))
    xvec = map(x -> x[2], filter(x -> x[1] > tstart, theta_t))

#    spl = fit(SmoothingSpline, tvec, xvec, lambda_spline)
#    xsmo = SmoothingSplines.predict(spl) # fitted vector
#    println("l=", length(xsmo))

#    if peakpar.sm && runpar.open_file_for_par
#        println(fl.wsk, "  ")
#        for ilength in 1:length(xsmo)
#            printfmtln(fl.wsk, "{1:f} {2:f}", tvec[ilength], xsmo[ilength])
#	end
#    end

    return tvec, xvec
end #smooth_whisking_trace(theta_t, peakpar, partohop, fl)

# This function finds the extrema of the whisking curve.
function find_extrema(tv, xv, whisk_set, whisk_amp, frac_sd_mnmx,
    list_of_extr, peakpar, partohop, awsval, fl)

    T_extr_one_cycle = 10.0

    global_min = findmin(xv)
    peakpar.sm ? println(fl.out, "global_min=", global_min, " ",
      tv[global_min[2]], " ", xv[global_min[2]]) : nothing
      
#   [t, x, i, +-1(max/min)] of the extremum.
    push!(list_of_extr, extr_val())

    list_of_extr[end].t = tv[global_min[2]]
    list_of_extr[end].x = global_min[1]
    list_of_extr[end].i = global_min[2]
    list_of_extr[end].mnmx = -1
    list_of_extr[end].t_extrapolate = -1.0

    println(fl.out, "min: ", list_of_extr[end])
     
    find_local_extrema(tv, xv, whisk_set, whisk_amp, list_of_extr, -1,
      frac_sd_mnmx, fl)
    reverse!(list_of_extr)
    find_local_extrema(tv, xv, whisk_set, whisk_amp, list_of_extr,  1,
      frac_sd_mnmx, fl)

#   Removing spurious extremum at the edges -  if the first maximum is smaller
#   than the value at the edge - it is not a real maximum.
#   println("tv1=", tv[1], " xv1=", xv[1])

    if list_of_extr[1].mnmx == 1 # maximum
        println("max remove item len=", length(list_of_extr))
        if list_of_extr[1].mnmx < xv[1]
	    for ilist in 1:length(list_of_extr)-1
	        list_of_extr[ilist] = deepcopy(list_of_extr[ilist+1])
	    end
	    deleteat!(list_of_extr, length(list_of_extr))
	    println("after remove len=", length(list_of_extr))
	end
    else
        list_of_extr[1].mnmx == -1 # minimum
        println("min remove item len=", length(list_of_extr))
        if list_of_extr[1].mnmx > xv[1]
	    for ilist in 1:length(list_of_extr)-1
	        list_of_extr[ilist] = deepcopy(list_of_extr[ilist+1])
	    end
	    deleteat!(list_of_extr, length(list_of_extr))
	    println("after remove len=", length(list_of_extr))
	end
    end
    
    print_extremum(list_of_extr, peakpar, fl)
    compute_average_maxima(list_of_extr, peakpar, awsval, fl)

end #find_extrema

# This function prints the variables of the extremum on the file fext.
function print_extremum(list_of_extr, peakpar, fl)

    if peakpar.sm && runpar.open_file_for_par
        println(fl.wsk, "  ")
        for iextr in 1:length(list_of_extr)
            extrval = list_of_extr[iextr]
            printfmtln(fl.wsk, "{1:f}  {2:f} {3:d} {4:d} {5:d}", extrval.t,
              extrval.x, extrval.i, extrval.mnmx, iextr)
        end
    end
    
end #print_extremum

# This function computes the average extremum.
function compute_average_maxima(list_of_extr, peakpar, awsval, fl)

    n_theta_max = 0
    theta_max = 0.0

#   Computing the maxima (set point)
    for iextr in 1:length(list_of_extr)
        extrval = list_of_extr[iextr]
#       printfmtln("{1:f}  {2:f} {3:d} {4:d} {5:d}", extrval.t,
#         extrval.x, extrval.i, extrval.mnmx, iextr)
	if extrval.mnmx == 1
	    n_theta_max += 1
	    theta_max += extrval.x
	end
    end

    if n_theta_max >= 1
        theta_max /= n_theta_max
    else
        theta_max = 0.0 #-999.9
    end

    printfmtln(fl.out, "n_theta_max={1:d} theta_max={2:f}", n_theta_max,
      theta_max)
#   printfmt(fl.aws, " {1:f}", theta_max)
    awsval.theta_max = theta_max

#   computing full amplitude
    iextr_beg = list_of_extr[1].mnmx == -1 ? 1 : 2
    iextr_end = list_of_extr[end].mnmx == 1 ? length(list_of_extr) - 1 :
      length(list_of_extr) - 2
    println("iextr_beg=", iextr_beg, " iextr_end=", iextr_end)
      
    n_theta_amp = 0
    theta_amp = 0.0
    if iextr_end >= iextr_beg
        for iextr in iextr_beg:2:iextr_end
	    extrval_min = list_of_extr[iextr]
	    extrval_max = list_of_extr[iextr+1]
	    n_theta_amp += 1
	    theta_amp += list_of_extr[iextr+1].x - list_of_extr[iextr].x
#	    println("theta_amp: ", n_theta_amp, " ", list_of_extr[iextr].x,
#	    " ", list_of_extr[iextr+1].x, " ", theta_amp)
	end

        if n_theta_amp >= 1
            theta_amp /= n_theta_amp
        else
            theta_amp = 0.0 #-999.9
	end
    else
        theta_amp = 0.0#-999.8
    end
    
    printfmtln(fl.out, "n_theta_amp={1:d} theta_amp={2:f}", n_theta_amp,
      theta_amp)
    awsval.theta_amp = theta_amp
    
end #compute_average_maxima

# This function finds the extrema of the smooth curves, on one side of the
# global minimum.
function find_local_extrema(tv, xv, x_av, x_sd, list_of_extr, direction,
         frac_sd_mnmx, fl)

    look_for = 1 #max

    beg_i_vec = min(list_of_extr[end].i + 2 , length(xv)-1)
    end_i_vec = direction==1 ? length(tv)-1 : 2
    println(fl.out, "end_i_vec=", end_i_vec)
#   println("beg_i_vec=", beg_i_vec, " end_i_vec=", end_i_vec, " direction=",
#     direction, " lxv=", length(xv), " ltv=", length(tv), " lle=",
#     length(list_of_extr))

    for iv in beg_i_vec : direction : end_i_vec
    
        extr_cond = look_for * (xv[iv] - xv[iv-1]) > 0 &&
	            look_for * (xv[iv] - xv[iv+1]) > 0 &&
		    look_for * (xv[iv] - list_of_extr[end].x) >
		    frac_sd_mnmx * x_sd

        extr_cond_opposite = -look_for * (xv[iv] - xv[iv-1]) > 0 &&
	                     -look_for * (xv[iv] - xv[iv+1]) > 0 &&
			     -look_for * (xv[iv] - list_of_extr[end].x) > 0

#      println(fl.out, iv, " ", extr_cond, " ", extr_cond_opposite, " ",
#        look_for, " ", tv[iv], " ", xv[iv-1], " ", xv[iv], " ", xv[iv+1],
#        " ", list_of_extr[end].x)
		    
        if extr_cond
	    push!(list_of_extr, extr_val())
            list_of_extr[end].t = tv[iv]
            list_of_extr[end].x = xv[iv]
            list_of_extr[end].i = iv
            list_of_extr[end].mnmx = -list_of_extr[end-1].mnmx
	    list_of_extr[end].t_extrapolate = -1.0
	    println(fl.out, "extr: ", list_of_extr[end])
	    look_for = -look_for 
	end

        if extr_cond_opposite
            list_of_extr[end].t = tv[iv]
            list_of_extr[end].x = xv[iv]
            list_of_extr[end].i = iv
	    list_of_extr[end].t_extrapolate = -1.0
	    println(fl.out, "ex_c: ", list_of_extr[end])
	end
    end

end #find_local_extrema

# This function computes the numbers and the start and end times of the
# whsiking cycles within each sniffing cycle.
function compute_whisking_cycles_in_sniffing_cycles(prebot_burst_start,
      list_of_extr, peakpar, partohop, fl)

#   i_prebot, n_peak_list_first, n_peak_list_last
#   whisking_within_breathing
    wskprop_ar = Array{wsk_properties,1}(undef,0)

#   println("prebot_burst_start=", prebot_burst_start)
#   println("list_of_extr=", list_of_extr)

    tstart = partohop.Tall - partohop.tstat - partohop.epsilon
    n_prebot_first_period = 1
    
    while n_prebot_first_period < length(prebot_burst_start) &&
      prebot_burst_start[n_prebot_first_period] < tstart
        n_prebot_first_period += 1
    end
    
    n_prebot_last_period = length(prebot_burst_start) - 1
    peakpar.sm ? println(fl.out, "n_prebot_first_period=",
      n_prebot_first_period, " n_prebot_last_period=", n_prebot_last_period) :
      nothing

#   Finding peak whisking within each breathing cycle.
#   This leaves one trough before the first peak, and the height can be
#   calculated.
    i_exter_first = list_of_extr[2].mnmx == 1 ? 2 : 3  #peak, not deep.
    
    for i_prebot in n_prebot_first_period:n_prebot_last_period
        while i_exter_first < length(list_of_extr)
            if list_of_extr[i_exter_first].t <= prebot_burst_start[i_prebot] &&
	      i_exter_first + 2 <= length(list_of_extr)
	        i_exter_first += 2
	    else
		break
	    end
        end
	
	i_exter_last = i_exter_first
	while i_exter_last+2 < length(list_of_extr) &&
	  list_of_extr[i_exter_last+2].t <= prebot_burst_start[i_prebot+1]
	    i_exter_last += 2
	end
        println("x len=", length(list_of_extr), " i_exter_last=", i_exter_last)
	
        wskprop = wsk_properties()
	wskprop.i_prebot = i_prebot
	wskprop.i_exter_first = i_exter_first
	wskprop.i_exter_last = i_exter_last
        wskprop.amp = Array{Float64,1}(undef,
	  div(i_exter_last-i_exter_first, 2) + 1)

        println("i_exter_first=", i_exter_first, " i_exter_last=",
	  i_exter_last, " len=", length(list_of_extr))
        for i_peak in i_exter_first:2:i_exter_last
	    amp = list_of_extr[i_peak].x - list_of_extr[i_peak-1].x
	    i_ar = 1 + div(i_peak-i_exter_first, 2)
	    wskprop.amp[i_ar] =  amp
        end
	
#       breathing_whisking_indices = [i_prebot, i_exter_first, i_exter_last]
	
	push!(wskprop_ar, wskprop)
#	whisking_within_breathing, breathing_whisking_indices)

        i_exter_first = i_exter_last
    end

    if peakpar.sm
        for wskprop in wskprop_ar
            printfmt(fl.out, "{1:d} {2:7.1f}", wskprop.i_prebot,
	      prebot_burst_start[wskprop.i_prebot])
	    printfmt(fl.out, " {1:d} {2:d} |", wskprop.i_exter_first,
	      wskprop.i_exter_last)
	    for i_amp in 1:length(wskprop.amp)
	        i_peak = wskprop.i_exter_first + 2 * (i_amp - 1)
	        printfmt(fl.out, " {1:7.1f}", list_of_extr[i_peak].t)
	    end
	    printfmt(fl.out, " |")
	    for i_amp in 1:length(wskprop.amp)
	        printfmt(fl.out, " {1:f}", wskprop.amp[i_amp])
	    end
	     printfmt(fl.out, "\n")
        end
#       println(fl.out, "wskprop_ar=", wskprop_ar)
    end

    i_peak_first = wskprop_ar[1].i_exter_first
    i_peak_last = wskprop_ar[end].i_exter_last
    i_peak_diff = i_peak_last - i_peak_first
    printfmtln(fl.out, "i_peak_first={1:d} i_peak_last={2:d} i_peak_diff={3:d}",
      i_peak_first, i_peak_last, i_peak_diff)

    if isodd(i_peak_diff)
        per_whisk = -999.1
    elseif i_peak_diff <= 0
        per_whisk = -999.0
    else
       T_whisk_interval = list_of_extr[i_peak_last].t -
          list_of_extr[i_peak_first].t
	T_per_whisk = T_whisk_interval * 2.0 / i_peak_diff
	printfmtln(fl.out, "T_whisk_interval={1:f} T_per_whisk={2:f}",
	  T_whisk_interval, T_per_whisk)
    end

    return wskprop_ar
end #compute_whisking_cycles_in_sniffing_cycles
    
# This function computes the statistics of whisking amplitudes within
# breathing cycles.
function statistics_of_wsk_in_breath(wskprop_ar, list_of_extr, peakpar,
    partohop, fl)
    
    npeak_ar = map(x -> length(x.amp), wskprop_ar)
    npeak_max = maximum(npeak_ar)

    n_amp = zeros(Int64, npeak_max)
    amp_av = zeros(npeak_max)
    amp_avt = zeros(npeak_max)
    amp_sd = zeros(npeak_max)
    time_av = zeros(npeak_max-1)
    time_avt = zeros(npeak_max-1)
    time_sd = zeros(npeak_max-1)

    for wskprop in wskprop_ar
        for i_amp in 1:length(wskprop.amp)
	    n_amp[i_amp] += 1
	    amp_av[i_amp] += wskprop.amp[i_amp]
	    amp_avt[i_amp] += wskprop.amp[i_amp] * wskprop.amp[i_amp]
	end
	
        for i_amp in 1:length(wskprop.amp)-1
	    i_peak_before = i_peak = wskprop.i_exter_first + 2 * (i_amp - 1)
	    t_before = list_of_extr[i_peak_before].t
	    i_peak_after = i_peak = wskprop.i_exter_first + 2 * (i_amp+1 - 1)
	    t_after = list_of_extr[i_peak_after].t
	    time_av[i_amp] += t_after - t_before
	    time_avt[i_amp] += (t_after - t_before)^2
	    printfmtln(fl.out, "i={1:d} t={2:f} {3:f} {4:f} {5:f}", i_amp,
	      t_before, t_after, t_after - t_before, time_av[i_amp])
	end
    end

    for i_amp in 1:npeak_max
        if n_amp[i_amp] > 1
	    amp_av[i_amp] /= n_amp[i_amp]
	    amp_avt[i_amp] /= n_amp[i_amp]
	    diff = amp_avt[i_amp] - amp_av[i_amp] * amp_av[i_amp]
	    if diff > 0.0
	        amp_sd[i_amp] = sqrt(diff)
	    elseif diff > -partohop.epsilon
	        amp_sd[i_amp] = 0.0
	    else
	        whisk_amp = -998.7
	    end
	elseif n_amp[i_amp] == 1	    
	    amp_sd[i_amp] = 0.0
	else
	    amp_av[i_amp] = -999.9
	    amp_sd[i_amp] = -999.8
	end
    end

    println(fl.out, "n_amp=", n_amp)
    println(fl.out, "amp_av=", amp_av)
    println(fl.out, "amp_sd=", amp_sd)

    for i_amp in 1:npeak_max-1
        if n_amp[i_amp] > 1
	   time_av[i_amp] /= n_amp[i_amp+1]
	   time_avt[i_amp] /= n_amp[i_amp+1]
	   diff = time_avt[i_amp] - (time_av[i_amp])^2
	   if diff > 0.0
	        time_sd[i_amp] = sqrt(diff)
	    elseif diff > -partohop.epsilon
	        time_sd[i_amp] = 0.0
	    else
	        time_amp = -988.7
	    end
	elseif n_amp[i_amp] == 1
	    time_sd[i_amp] = 0.0
	else
	    time_av[i_amp] = -989.9
	    time_sd[i_amp] = -989.8
	end
    end
    
    println(fl.out, "time_av=", time_av)
    println(fl.out, "time_sd=", time_sd)

# Statistics of timing

    return amp_av, time_av
end #statistics_of_wsk_in_breath

# This function computes the time from the breathing onset to first whisking
# peak afterwards vs. the time from the sniffing onset to the last whisking
# peak of the previous breathing cycle.
# The function also computes the linear regression of the T_diff_post vs.
# T_diff_pre curve, for T_diff_pre in [T_diff_pre_min, T_diff_pre_max].
function timing_of_last_first_whisking(prebot_burst_start, wskprop_ar,
    list_of_extr, peakpar, partohop, awsval, fl)

    T_diff_pre_min = runpar.T_diff_pre_interval[1] #40.0
    T_diff_pre_max = runpar.T_diff_pre_interval[2] #120.0
    
    n_T_diff_inside_range = 0
    T_diff_pre_ar = Array{Float64,1}(undef,0)
    T_diff_whisk_ar = Array{Float64,1}(undef,0)

    println(fl.out, "\nT_diff  len=", length(wskprop_ar))
    amp_avr  = 0.0
    amp_avrt = 0.0
    
    for iwskprop in 3:length(wskprop_ar)
        wskprop_early = wskprop_ar[iwskprop-1]
        wskprop_now = wskprop_ar[iwskprop]
	T_breath = prebot_burst_start[wskprop_now.i_prebot]

	i_earlier_peak = wskprop_early.i_exter_first +
	  2 * (length(wskprop_early.amp) - 1)
	T_earlier_whisk = list_of_extr[i_earlier_peak].t

        i_later_peak = wskprop_now.i_exter_first
	T_later_whisk = list_of_extr[i_later_peak].t
	
        T_diff_pre = T_breath - T_earlier_whisk
	T_diff_post = T_later_whisk - T_breath

        T_diff_whisk = T_later_whisk - T_earlier_whisk

        amp_avr  += wskprop_now.amp[1]
        amp_avrt += wskprop_now.amp[1] * wskprop_now.amp[1]
	
        if peakpar.sm && runpar.open_file_for_par
            printfmt(fl.tpp, "{1:f} {2:f} {3:f} {4:d} {5:f} {6:f} {7:f}",
	      T_diff_pre, T_diff_post, T_diff_whisk,
	      wskprop_now.i_prebot, T_earlier_whisk, T_breath, T_later_whisk)
	    printfmtln(fl.tpp, " {1:f}", wskprop_now.amp[1])
        end
	
        if T_diff_pre > T_diff_pre_min && T_diff_pre < T_diff_pre_max
	    n_T_diff_inside_range += 1
	    push!(T_diff_pre_ar, T_diff_pre)
	    push!(T_diff_whisk_ar, T_diff_whisk)
	end	
    end

    if length(wskprop_ar) >= 2
        amp_avr  /= length(wskprop_ar) - 1
        amp_avrt /= length(wskprop_ar) - 1
	diff = amp_avrt - amp_avr * amp_avr
	if diff > 0.0
	    amp_sd = sqrt(diff)
	elseif diff > -partohop.epsilon
	    amp_sd = 0.0
	else
	    amp_sd = -999.7
	end
    else
        amp_avr = -999.9
	amp_sd  = -999.8
    end	

#   computing linear regression
    Mmat = zeros(length(T_diff_pre_ar),2)
    Mmat[:,1] = T_diff_pre_ar
    Mmat[:,2] = ones(length(T_diff_pre_ar))
    coeff_pred = Mmat \ T_diff_whisk_ar
    println(fl.out, "coeff_pred=", coeff_pred)
    
#   lin_fit_Tdiff = Polynomials.fit(T_diff_pre_ar, T_diff_whisk_ar, 1)
    printfmtln(fl.out, "n_T_diff_inside_range={1:d}", n_T_diff_inside_range)
    
#   printfmtln(fl.out, "Linear regression coefficients: {1:f} {2:f}",
#     lin_fit_Tdiff[0], lin_fit_Tdiff[1])
#   println(fl.out, " lin_fit_Tdiff=", lin_fit_Tdiff)
#   printfmt(fl.aws, " {1:f} {2:f}", lin_fit_Tdiff[0], lin_fit_Tdiff[1])
    printfmtln(fl.out, "Linear regression coefficients: {1:f} {2:f}",
      coeff_pred[2], coeff_pred[1])
#   printfmt(fl.aws, " {1:f} {2:f} {3:f} {4:f}", coeff_pred[2], coeff_pred[1],
#     amp_avr, amp_sd)
    awsval.coeff_pred2 = coeff_pred[2]
    awsval.coeff_pred1 = coeff_pred[1]
    awsval.amp_avr = amp_avr
    awsval.amp_sd = amp_sd
    
end #timing_of_last_first_whisking

# This function computs whisking amplitudes and the statistics of whisking
# angles within a breathing cycle.
function whisking_within_sniffing(theta_t, whisk_set, whisk_amp, peakpar,
    prebot_burst_start, partohop, awsval, fl)

#   prebot_burst_start = Array{Float64,1}(undef,0)
    prot_in_sniff_start = Array{Array{Any,1}}(undef,0)

    list_of_extr = Array{extr_val,1}(undef,0)

#   println("rstpbstr.prebot_burst_start=", rstpbstr.prebot_burst_start)
#   find_time_external_prebot(prebot_burst_start, peakpar, partohop, fl)

    tvec, xvec = extract_smooth_whisking_trace(theta_t, peakpar, partohop, fl)

    find_extrema(tvec, xvec, whisk_set, whisk_amp, partohop.frac_sd_mnmx,
      list_of_extr, peakpar, partohop, awsval, fl)

    if length(list_of_extr) > 3
        wskprop_ar = compute_whisking_cycles_in_sniffing_cycles(
	  prebot_burst_start, list_of_extr, peakpar, partohop, fl)

        height_maxima, del_time_maxima = statistics_of_wsk_in_breath(wskprop_ar,
	  list_of_extr, peakpar, partohop, fl)

        timing_of_last_first_whisking(prebot_burst_start, wskprop_ar,
	  list_of_extr, peakpar, partohop, awsval, fl)
    else
        height_maxima = -999.9
	del_time_maxima = -998.9
    end
    
    return height_maxima, del_time_maxima
end #whisking_within_sniffing

# This function prints the values within awsval on the file fl.aws
function aws_print(awsval, fl)

    printfmt(fl.aws, "{1:f} {2:f}", awsval.whisk_set, awsval.whisk_amp)
    
    if abs(awsval.coeff_pred2) < 1.0e10
        printfmt(fl.aws, " {1:f}", awsval.coeff_pred2)
    else
        printfmt(fl.aws, " {1:f}", -999.989)
    end
    
    if abs(awsval.coeff_pred1) < 1.0e10
        printfmt(fl.aws, " {1:f}", awsval.coeff_pred1)
    else
        printfmt(fl.aws, " {1:f}", -999.989)
    end
    
    if abs(awsval.amp_avr) < 1.0e10
        printfmt(fl.aws, " {1:f}", awsval.amp_avr)
    else
        printfmt(fl.aws, " {1:f}", -999.989)
    end
    
    if abs(awsval.amp_sd) < 1.0e10
        printfmt(fl.aws, " {1:f}", awsval.amp_sd)
    else
        printfmt(fl.aws, " {1:f}", -999.989)
    end
    printfmt(fl.aws, " {1:f} {2:f}", awsval.theta_max, awsval.theta_amp)

end #aws_print

# This function computes the consecutive whsiking amplitudes.
function find_whisking_amplitudes(rstpbstr, sval, partohop, suffix, mode, fl)
    
    peakpar = peak_par()
    peakpar.deltat_spk = 0.05
    peakpar.deltat = 0.1
    peakpar.nburst = 3  # Does not consider the first burst.
    peakpar.sm = true #false #true

    awsval = aws_val()

    theta_t = Array{Array{Float64,1}}(undef,0)

#   Computing the whisking angle based on motorneuron firing.
    whisk_set, whisk_amp = compute_write_whisking_angle(theta_t, peakpar,
      rstpbstr.rst_ar, partohop, awsval, fl)

#   computing the statistics of whisking angles within a breathing cycle.
    height_maxima, del_time_maxima = whisking_within_sniffing(theta_t,
      whisk_set, whisk_amp, peakpar, rstpbstr.prebot_burst_start, partohop,
      awsval, fl)

    aws_print(awsval, fl)
    printfmt(fl.aws, "\n")
    flush(fl.aws)

    return height_maxima, del_time_maxima
end
