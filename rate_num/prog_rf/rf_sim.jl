# This file contains the functions for simulating the dynamical system.

# This function simulates a network with one parameter set.
function one_par(suffix, netpar, sval, av, fl)

    pbstr = pb_str()
    pbstr.s_preBot = 0.0
    pbstr.s_preBot_o = 0.0
    pbstr.up = false
    pbstr.n_preBot_cycle = 0
    pbstr.m_dat_cycle = 3 # tmin, tmax, T_up
    m_TB = ceil(Int64, runpar.Tall / netpar.TB)
    println(fl.out, "n_TB=", m_TB)
    pbstr.dat_cycle = Array{Array{Float64,1},1}(undef, m_TB)
#   Tmin, Tmax, Mr, Mf    
    [pbstr.dat_cycle[i_TB] = zeros(4) for i_TB in 1:m_TB]

    statstr = stat_str()

    substitute_values_in_structures(netpar, runpar, statstr, av, fl)
    substitute_parameters_in_netpar(netpar, runpar, sval, fl)

    Varbar = Array{Float64,1}(undef, 2*netpar.Npop+1)

    #Initial conditions
    in_con(Varbar, netpar, runpar, fl)

    store_Varbar = Array{Array{Float64,1},1}(undef, 0) #2*netpar.Npop+1)

    #Integrating the differential euqations over the time interval
    n_run(Varbar, pbstr, netpar, runpar, statstr, store_Varbar, fl)

    #Computing statistics of activity
    statistics_of_activity(statstr, pbstr, netpar, runpar, av, fl)

    #Computing statistics of whisking
    statistics_of_whisking(statstr, pbstr, netpar, runpar, av, fl)

    #Finding the temporal values for starting and end of theta values
    netpar.grB > runpar.epsilon ? processed_V_values(store_Varbar, fl) : nothing

end #one_par()

# This function substitues values into the structures netpar and runpar.
function substitute_values_in_structures(netpar, runpar, statstr, av, fl)

    runpar.NT = round(Int64, runpar.Tall/runpar.deltat)
    runpar.tstat > runpar.Tall ? runpar.tstat = runpar.Tall : nothing

    if runpar.smforce == 'p'
        runpar.sm = true
    elseif runpar.smforce == 'n'
        runpar.sm = false
    elseif runpar.smforce == 'l'
        if sval.scan_type == 'n'
	    runpar.sm = true
	else
	    runpar.sm = false
	end
    else
        println("smforce should be either p or n or l or a !!! smforce=",
        runpar.smforce)
	exit(0)
    end
    println(fl.out, "sm=", runpar.sm)
    flush(fl.out)

    statstr.hist_store_ar = zeros(runpar.n_hist_store)
    statstr.peak_time_ar  = Array{Float64,1}(undef,0)
    statstr.Mav  = zeros(netpar.Npop)
    statstr.Mavt = zeros(netpar.Npop)
    statstr.Msd  = zeros(netpar.Npop)

    statstr.theta_cb = CircularBuffer{Float64}(runpar.n_hist_store)
    fill!(statstr.theta_cb, 0)
    statstr.theta_min_ar = Array{Array{Float64,1},1}(undef,0)
    statstr.theta_max_ar = Array{Array{Float64,1},1}(undef,0)

    statstr.theta_sum = 0.0
    statstr.theta_sumt = 0.0
    
    av.Mav  = zeros(netpar.Npop)
    av.Msd  = zeros(netpar.Npop)
        
    av.Mup  = zeros(netpar.Npop)

end #substitute_values_in_structures

# This function substitues parameter into netpar.
function substitute_parameters_in_netpar(netpar, runpar, sval, fl)

    netpar.Ja = 1000.0 * netpar.gamma * netpar.gKm 
    netpar.Jintra  = netpar.gintra * netpar.DelV / netpar.taus
    netpar.Jinter = netpar.ginter * netpar.DelV / netpar.taus

    printfmtln(fl.out, "Ja={1:f} Jintra={2:f} Jinter={3:f}", netpar.Ja,
      netpar.Jintra, netpar.Jinter)

    netpar.IrB = netpar.grB * netpar.DelV
    printfmtln(fl.out, "IrB={1:f}", netpar.IrB)

    if sval.scan_type == 'n'
        printfmtln(fl.grb, "-20.0 0.0")
        printfmtln(fl.grb, "0.0 0.0")
        printfmtln(fl.grb, "0.0 {1:f}", netpar.grB)
        printfmtln(fl.grb, "{1:f} {2:f}", netpar.Tup, netpar.grB)
        printfmtln(fl.grb, "{1:f} 0.0", netpar.Tup)
        printfmtln(fl.grb, "120.0 0.0")
    end
    
    
    netpar.Ja_F = 1000.0 * netpar.gamma_F * netpar.gKm_F 
    netpar.JFr = netpar.gFr * netpar.DelV_F / netpar.taus

    printfmtln(fl.out, "Ja_F={1:f} JFr={2:f}", netpar.Ja_F, netpar.JFr)

end #substitute_parameters_in_netpar

#This function substitutes the initial conditions.
function in_con(Varbar, netpar, runpar, fl)

    if (runpar.incond == 'r') # read initial conditions
        for ipop in 1:netpar.Npop, ivar in 1:2
	    Varbar[2*(ipop-1) + ivar] = netpar.initial_val[2*(ipop-1) + ivar]
	end
	Varbar[2*netpar.Npop + 1] = netpar.initial_val[2*netpar.Npop + 1]
    else
        println("wrong incond=", runpar.incond)
	exit(0)
    end
    
end #in_con

# This function solves the ODEs.
function n_run(Varbar, pbstr, netpar, runpar, statstr, store_Varbar, fl)

    nvar = 2*netpar.Npop+1 # Npop of [s,a] + theta

    k0 = zeros(nvar)
    k1 = zeros(nvar)
    k2 = zeros(nvar)
    k3 = zeros(nvar)
    k4 = zeros(nvar)
    Varc = zeros(nvar)
    
    it = 0
    
    time = it * runpar.deltat
    if runpar.tmcol + runpar.epsilon >= runpar.Tall
        pr_fct(Varbar, pbstr, netpar, runpar, time, it, store_Varbar, fl)
    end

    for it = 1:runpar.NT
        time = it * runpar.deltat

        compute_s_preBot(pbstr, time, it, netpar, runpar, fl)

        #Advancing the state vetor for one time step.
	if runpar.method  == 'r'                #Runge-Kutta-4 method
            one_integration_step(Varbar, k0, k1, Varc, 0.0            ,
	      time, it, pbstr, netpar, runpar, fl)
            one_integration_step(Varbar, k1, k2, Varc, runpar.deltat/2,
	      time, it, pbstr, netpar, runpar, fl)
            one_integration_step(Varbar, k2, k3, Varc, runpar.deltat/2,
	      time, it, pbstr, netpar, runpar, fl)
            one_integration_step(Varbar, k3, k4, Varc, runpar.deltat  ,
	      time, it, pbstr, netpar, runpar, fl)
	    Varbar += (runpar.deltat/6.0) * (k1 + 2.0 * (k2 + k3) + k4)
	    
        elseif runpar.method  == 't'            #Runge-Kutta-2 method
            one_integration_step(Varbar, k0, k1, Varc, 0.0            ,
	      time, it, pbstr, netpar.Ntot, netpar, runpar, fl)
            one_integration_step(Varbar, k1, k2, Varc, runpar.deltat/2,
	      time, it, pbstr, netpar, runpar, fl)
            Varbar += runpar.deltat * k2

        elseif runpar.method  == 'e'            #Euler method 
            one_integration_step(Varbar, k0, k1, Varc, 0.0            ,
	      time, it, pbstr, netpar, runpar, fl)
            Varbar += runpar.deltat * k1
        end

#       Checking for nan.
	for ivar in 1:nvar
	    if any(isnan, Varbar[ivar])
	        println("NaN in Varbar[ivar]!", var)
	        exit(0)
	    end
	end

	if time >= runpar.Tall -runpar.tmcol + runpar.epsilon &&
	  it%runpar.twrite == 0
            pr_fct(Varbar, pbstr, netpar, runpar, time, it, store_Varbar, fl)
	end
	
	store_history(Varbar, statstr, pbstr, netpar, runpar, time, it, fl)
	if it > runpar.n_hist_store &&
	  time > runpar.Tall - runpar.tstat - runpar.epsilon
	    find_maxima(statstr, netpar, runpar, time, it, fl)
	end
	
#       Storing theta in circular buffer
        push!(statstr.theta_cb, Varbar[2*netpar.Npop+1])
#       Finding extrema in theta
        if time >= runpar.Tall -runpar.tstat + runpar.epsilon &&
	  it >= runpar.n_hist_store
	    find_theta_exremum(time, netpar, runpar, statstr, fl)
        end	
    end

#   Normalizing M values during preBot activity.
#   for i_TB in 1:length(pbstr.dat_cycle)
#       [pbstr.dat_cycle[i_TB][iM] *= 1000.0 / netpar.Tup for iM in 3:4]
#   end

#   println("pbstr.dat_cycle=", pbstr.dat_cycle)
end #n_run

# This function computes s_preBot based on the time and the preBot phase.
function compute_s_preBot(pbstr, time, it, netpar, runpar, fl)

    pbstr.s_preBot_o = pbstr.s_preBot
    pbstr.Tphase = mod(time - runpar.epsilon, netpar.TB)
    pbstr.s_preBot = pbstr.Tphase < netpar.Tup ? 1.0 : 0.0
#   printfmtln("so={1:e} s={2:e} it={3:d} time={4:f}", pbstr.s_preBot_o,
#     pbstr.s_preBot, it, time)

    if pbstr.s_preBot_o <= runpar.epsilon && pbstr.s_preBot > runpar.epsilon
        pbstr.up = true
	pbstr.n_preBot_cycle += 1
	if pbstr.n_preBot_cycle > length(pbstr.dat_cycle)
	    println("n_preBot_cycle=", pbstr.n_preBot_cycle, " > ",
	      length(pbstr.dat_cycle))
	    exit(0)
	end
#       println(fl.out, "it=", it, " time=", time, " up n=",
#	  pbstr.n_preBot_cycle)
#       Tmin	  
	pbstr.dat_cycle[pbstr.n_preBot_cycle][1] = time - runpar.deltat
    elseif pbstr.s_preBot_o > runpar.epsilon && pbstr.s_preBot <= runpar.epsilon
        pbstr.up = false
#       println(fl.out, "it=", it, " time=", time, " down")
#       Tmax
	pbstr.dat_cycle[pbstr.n_preBot_cycle][2] = time - runpar.deltat
    end

end #compute_s_preBot

# This function computes M.
function Mcal(ss, aa, time, pbstr, netpar, runpar, fl)
#   Linear-threshold function	 
    linthr(xx) = xx > 0 ? xx : 0
#   bb = 100  
#   linthr(x)=x/(1.0+exp(-bb*x))

#   Tphase = mod(time - runpar.epsilon, netpar.TB)
#   s_preBot = Tphase < netpar.Tup ? 1.0 : 0.0
    I_preBot_to_vIRtr = netpar.IrB * pbstr.s_preBot

    MM = zeros(netpar.Npop)
    
    for ipop in 1:2   # vIRt_r, vIRt_p
        inhibition_factor = ipop == 1 ? netpar.inhibition_factor : 1.0
        I_input = netpar.I0 - netpar.I0_threshold - aa[ipop] -
	  inhibition_factor * (netpar.Jintra * ss[ipop] + netpar.Jinter * 
	  ss[3-ipop])
#       preBot->vIRtr inhibition.
	I_input -= ipop == 1 ? I_preBot_to_vIRtr * inhibition_factor : 0.0
        MM[ipop] = netpar.beta * linthr(I_input)
    end

    ipop = netpar.Npop   # FN
    I_input_F = netpar.I0_F - netpar.I0_threshold_F - aa[ipop] -
      netpar.JFr * ss[1]
    MM[ipop] = netpar.beta_F * linthr(I_input_F)

#   println("ss=", ss, " aa=", aa, " MM=", MM)

    return MM
end #Mcal


# This function computes the force term from MF, using the approximation
# computed in the program cfit.jl
function force_from_MF(MF, netpar, runpar, fl)

    AA = 9.225334
          
    f1 = 525.896192
    B2 = -23.013849
    f2 = 611.851637
    
    B3 = 152.398764
    f3 = 459.910247
    
    AL = 1.020052
    freql = 76.996960

    freq = 1000.0 * MF
    zz = 1000.0 * MF / freql

    xlog = AL * log(1.0 + zz)
    poly_freq = (freq / f1) + B2 * (freq / f2)^2 + B3 * (freq / f3)^3
    fit1 = AA * poly_freq / (1 + poly_freq)
    fit_val = fit1 + xlog

    return fit_val
end #force_from_MF

#This function computes one integration step.
function one_integration_step(Varbar, kin, kout, Varc, delt, time, it,
         pbstr, netpar, runpar, fl)

#   Runge-Kutta input variables
    Varc = delt > runpar.epsilon ? Varc = Varbar + delt * kin : Varbar
    MM = zeros(netpar.Npop)
    
    ss = [Varc[(ipop-1)*2+1] for ipop in 1:netpar.Npop]
    aa = [Varc[(ipop-1)*2+2] for ipop in 1:netpar.Npop]
    MM = Mcal(ss, aa, time, pbstr, netpar, runpar, fl)
    theta = Varc[2*netpar.Npop+1]
    
    for ipop in 1:2   # vIRt_r, vIRt_p
# s
        kout[(ipop-1)*2+1] = -ss[ipop] / netpar.taus + MM[ipop]
# a
        kout[(ipop-1)*2+2] = (-aa[ipop] + netpar.Ja * MM[ipop]) /
	  netpar.taua
    end

    ipop = netpar.Npop   # FN
# s
    kout[(ipop-1)*2+1] = -ss[ipop] / netpar.taus + MM[ipop]
# a
    kout[(ipop-1)*2+2] = (-aa[ipop] + netpar.Ja_F * MM[ipop]) /
	  netpar.taua_F

    force_term = force_from_MF(MM[3], netpar, runpar, fl)
#   printfmtln(fl.out, "time={1:f} MF={2:f} force_term={3:f}", time, MM[3],
#     force_term)
    kout[2*netpar.Npop+1] = -(theta / netpar.tauw) + force_term
##  kout[2*netpar.Npop+1] = -(theta / netpar.tauw) + 50.0 * MM[3]

end #one_integration_step()

#This function writes the data as functions of time
function pr_fct(Varbar, pbstr, netpar, runpar, time, it, store_Varbar, fl)

    ss = [Varbar[(ipop-1)*2+1] for ipop in 1:netpar.Npop]
    aa = [Varbar[(ipop-1)*2+2] for ipop in 1:netpar.Npop]
    MM = Mcal(ss, aa, time, pbstr, netpar, runpar, fl)
    
    if runpar.sm
        printfmt(fl.col, "{1:f}", time)
        for ipop = 1:netpar.Npop
           printfmt(fl.col, " {1:11.8f} {2:f} {3:f}", MM[ipop], ss[ipop],
	     aa[ipop])
        end
        printfmt(fl.col, " {1:f}", Varbar[2*netpar.Npop+1])
        printfmt(fl.col, "\n")
    end

    push!(store_Varbar, [time, Varbar[2*netpar.Npop+1]])
        
end #pr_fct

# This function stores the history of variables.
function store_history(Varbar, statstr, pbstr, netpar, runpar, time, it, fl)

    ss = [Varbar[(ipop-1)*2+1] for ipop in 1:netpar.Npop]
    aa = [Varbar[(ipop-1)*2+2] for ipop in 1:netpar.Npop]
    MM = Mcal(ss, aa, time, pbstr, netpar, runpar, fl)

    i_var_store = 1
    xx = MM[i_var_store]
    deleteat!(statstr.hist_store_ar, 1)
    push!(statstr.hist_store_ar, xx)

    if time > runpar.Tall - runpar.tstat - runpar.epsilon
        statstr.Mav  += runpar.deltat * MM
        statstr.Mavt += runpar.deltat * MM .* MM
#       println("it=", it, " MM=", MM, " Mav=", statstr.Mav, " Mav2=",
#         statstr.Mavt)
    end

    if pbstr.up
#       vIRt_r
        pbstr.dat_cycle[pbstr.n_preBot_cycle][3] += MM[1] * runpar.deltat
#       FN
        pbstr.dat_cycle[pbstr.n_preBot_cycle][4] += MM[3] * runpar.deltat
    end 

    if time > runpar.Tall - runpar.tstat + runpar.epsilon
        theta = Varbar[2*netpar.Npop+1]
        statstr.theta_sum += theta * runpar.deltat
	statstr.theta_sumt += theta * theta * runpar.deltat
#	println(fl.out, "theta=", theta, " theta_sum=", statstr.theta_sum, " theta_sumt=", statstr.theta_sumt)
    end
    
end #store_history

# This function finds the maxima of a variable over time.
function find_maxima(statstr, netpar, runpar, time, it, fl)
	      
    ndiv = div(runpar.n_hist_store, 2)
    i_argmax = argmax(statstr.hist_store_ar)
    if i_argmax == ndiv
#       time_peak = time - ndiv * runpar.deltat
        y0 = statstr.hist_store_ar[ndiv+1]
	y1 = statstr.hist_store_ar[ndiv]
	y2 = statstr.hist_store_ar[ndiv-1]
        time_peak, y_peak = detect_peak(y0, y1, y2, time, runpar.deltat,
	  runpar, fl)
#       printfmtln(fl.out, "{1:f} {2:f} {3:d}", time_peak,
#         statstr.hist_store_ar[i_argmax], it - ndiv)
        push!(statstr.peak_time_ar, time_peak)
    end

end #find_maxima

# This function computes the peak time.
function detect_peak(y0, y1, y2, time, deltat, runpar, fl)

    ndiv = div(runpar.n_hist_store, 2)
    
    if y1 >= y0 && y1 >= y2
        xb = y2 - y0
        xc = y0 - 2.0 * y1 + y2
	if (abs(xc) < runpar.epsilon)
            tpeak = time - runpar.deltat
	    ypeak = y1
        else
	    tpeak = time - (ndiv+1) * runpar.deltat +
	      0.5 * (xb / xc) * runpar.deltat
	    ypeak = y1 - 0.125 * xb * xb / xc
#	    printfmtln("y0={1:11.8f} y1={2:11.8f} y2={3:11.8f} tpeak={4:11.8f}"*
#	      " ypeak={5:f}", y0, y1, y2, tpeak, ypeak)
	end
    else
        println("no peak! y0=", y0, " y1=", y1, " y2=", y2)
	exit(0)
    end

    return tpeak, ypeak
end #detect_peak

# This function computes the statistics of firing activity.
function statistics_of_activity(statstr, pbstr, netpar, runpar, av, fl)

#   Finding if there is a fixed point.

    statstr.Mav  /= runpar.tstat + runpar.deltat
    statstr.Mavt /= runpar.tstat + runpar.deltat

#   println("tt=", runpar.tstat + 1, " Mav=", statstr.Mav, " Mav2=",
#     statstr.Mavt)
    for ipop in 1:netpar.Npop
        diff = statstr.Mavt[ipop] - (statstr.Mav[ipop])^2
	if diff > 0.0
	    statstr.Msd[ipop] = sqrt(diff)
	elseif diff > -runpar.epsilon
	    statstr.Msd[ipop] = 0.0
	else
	    println("ipop=", ipop, " diff=", diff, " < 0!")
	    exit(0)
	end
    end
    println(fl.out, "Msd=", statstr.Msd)

    av.Mav = deepcopy(statstr.Mav)
    av.Msd = deepcopy(statstr.Msd)

    epsilon_sd = 1.0e-6
#   Searching for a fixed point
    if av.Msd[1] < epsilon_sd && av.Msd[2] < epsilon_sd
        if av.Mav[1] > epsilon_sd && av.Mav[2] > epsilon_sd
            Tper = -998.3
	elseif av.Mav[1] > epsilon_sd
	    Tper = -998.1
	elseif av.Mav[2] > epsilon_sd
	    Tper = -998.2
	else
	    Tper = -998.0
	end
    else
##      println(fl.out, "list of peaks")
##      [printfmtln(fl.out, "{1:f}", tpeak) for tpeak in statstr.peak_time_ar]

        if length(statstr.peak_time_ar) >= 2
            Tper = (statstr.peak_time_ar[end] - statstr.peak_time_ar[1]) /
              (length(statstr.peak_time_ar) - 1)
	else
	    Tper = -999.5
	end
    end
    printfmtln(fl.out, "Tper={1:f}", Tper)
    av.Tper = Tper

    for i_TB in 3:length(pbstr.dat_cycle)
        av.Mup[1] += pbstr.dat_cycle[i_TB][3]
        av.Mup[3] += pbstr.dat_cycle[i_TB][4]
    end

    if (length(pbstr.dat_cycle) >= 3)
        av.Mup[1] /= length(pbstr.dat_cycle) - 2
        av.Mup[3] /= length(pbstr.dat_cycle) - 2
   end
   
end #statistics_of_activity

# This function computes the statistics of firing activity.
function statistics_of_whisking(statstr, pbstr, netpar, runpar, av, fl)

    compute_theta_av_sd(statstr, pbstr, netpar, runpar, av, fl)

    if av.theta_sd > runpar.min_theta_for_osci
        compute_whisk_amp_avr(statstr, pbstr, netpar, runpar, av, fl)
    else
        av.whisk_amp_avr = 0.0
        av.whisk_setup = av.theta_av
    end

end #statistics_of_whisking

# This function computes the average and the standard deviation of theta.
function compute_theta_av_sd(statstr, pbstr, netpar, runpar, av, fl)

    statstr.theta_sum /= runpar.tstat
    statstr.theta_sumt /= runpar.tstat
#   println(fl.out, "tstat=", runpar.tstat, " theta_sum=", statstr.theta_sum,
#     " theta_sumt=", statstr.theta_sumt)
    diff = statstr.theta_sumt - statstr.theta_sum^2
    theta_sd = 0.0
    if diff > 0.0
        theta_sd = sqrt(diff)
    elseif diff > -runpar.epsilon
        theta_sd = 0.0
    else
        theta_sd = -999.2
    end
    av.theta_av = statstr.theta_sum
    av.theta_sd = theta_sd
    
end #compute_theta_av_sd

# This function comutes the amplitude and the set point of theta.
function compute_whisk_amp_avr(statstr, pbstr, netpar, runpar, av, fl)

    if length(statstr.theta_min_ar) >=1 && length(statstr.theta_max_ar) >=1
        if runpar.sm
            for extremum_ar in (statstr.theta_min_ar, statstr.theta_max_ar)
                if length(extremum_ar) > 0
                    for (index, record) in enumerate(extremum_ar)
	                printfmtln(fl.pek, "{1:f} {2:f} {3:d}", record[1],
		        record[2], index)
	            end
	            printfmtln(fl.pek, "  ")
	        end
            end
        end

        i_beg_min = 1
        i_beg_max = statstr.theta_max_ar[1][1] > statstr.theta_min_ar[1][1] ?
	  1 : 2
        if statstr.theta_max_ar[i_beg_max][1] <
	  statstr.theta_max_ar[i_beg_min][1]
            printfmtln("t_max_2={1:f} > t_min_1={2:f} !")
	    exit(0)
        end
        n_whisks = length(statstr.theta_max_ar) + 1 - i_beg_max
        println(fl.out, "i_beg_min=", i_beg_min, " i_beg_max=", i_beg_max,
          " n_whisks=", n_whisks)
    else
        n_whisks = 0
    end

    whisk_amp_avr = 0.0
    whisk_setup = 0.0
    
    if n_whisks >= 1
        println(fl.out, "\nwhisk_amp")
    
        for iwhisk in i_beg_min:n_whisks
            t_min = statstr.theta_min_ar[iwhisk][1]
            theta_min = statstr.theta_min_ar[iwhisk][2]
	    iwhisk_max = iwhisk + i_beg_max - 1
 	    t_max = statstr.theta_min_ar[iwhisk][1]
            theta_max = statstr.theta_max_ar[iwhisk_max][2]

            if t_min > t_max
	        println("t_min=", t_min, " > t_max=", t_max, "!")
	        exit(0)
	    end
	
	    whisk_amp = theta_max - theta_min
            whisk_amp_avr += whisk_amp
            whisk_setup += theta_max
	    
#	    printfmtln(fl.out, "iwhisk={1:d} iwhisk_max={2:d} theta_min={3:f}" *
#	      " theta_max={4:f} whisk_amp={5:f} whisk_amp={6:f}", iwhisk,
#	      iwhisk_max, theta_min, theta_max, whisk_amp, whisk_amp_avr)
        end

#       println(fl.out, "nm=", n_whisks - i_beg_min + 1)
        whisk_amp_avr /= n_whisks - i_beg_min + 1
	whisk_setup /= n_whisks - i_beg_min + 1
    end
    printfmtln(fl.out, "whisk_amp_avr={1:f} whisk_setup={2:f}", whisk_amp_avr,
       whisk_setup)

    av.whisk_amp_avr = whisk_amp_avr
    av.whisk_setup = whisk_setup
    
end #compute_whisk_amp_avr

# This function find extrama in theta.
function find_theta_exremum(time, netpar, runpar, statstr, fl)

    ndiv = div(runpar.n_hist_store, 2)

#   Checking for a maximum
    check_for_extremum(statstr.theta_cb, 1, ndiv, statstr.theta_max_ar, time,
      runpar.deltat, runpar, fl)

#   Checking for a minimum
    check_for_extremum(-statstr.theta_cb, -1, ndiv, statstr.theta_min_ar, time,
      runpar.deltat, runpar, fl)
      
end #find_theta_exremum

# This function checks for an extremum
function check_for_extremum(array, mnmx, ndiv, str_ar_to_push, time, deltat,
  runpar, fl)

    i_argmax = argmax(array)
    if i_argmax == ndiv
        y0 = array[ndiv+1]
	y1 = array[ndiv]
	y2 = array[ndiv-1]
        time_peak, y_peak = detect_peak(y0, y1, y2, time, deltat, runpar, fl)
        push!(str_ar_to_push, [time_peak, mnmx * y_peak])
    end
    
end #check_for_extremum

# This function defines the sigmoid function Gammaf
function Gammaf(VV, theta, sigma)
    return 1.0/(1.0+exp(-(VV-theta)/sigma))
end #Gammaf
