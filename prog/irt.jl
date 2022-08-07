# This file contains the functions for running the simulation program.

#@time 
import Random

using Formatting
using Printf
using SmoothingSplines
using Statistics
using Polynomials
using DSP
using NLsolve
using StatsBase
#using Random

# This function sets the value of the variable sval.parname to the value
# sval.par.
function modify_parameter(svalx)

  str_command = svalx.parname[1] * " = " * string(svalx.par)
  println("str_command=", str_command)
  eval(Meta.parse(str_command))

  for str_command in svalx.parname[2:end]
      println("str_command=", str_command)
      eval(Meta.parse(str_command))
  end
  
end #modify_parameter

# This function opens files and names them according to scan_type.
function open_files(scan_type, suffix, ipara, ipar, irepeat, fl)

    if scan_type == 'n'
        name_end = ""
    elseif scan_type == 'y'
        name_end = "." * string(ipar) * "." * string(irepeat)
    elseif scan_type == 't'
         name_end = "." * string(ipara) * "." * string(ipar) * "." *
          string(irepeat)
    else
        println("wrong scan_type=", scan_type)
	exit(0)
    end
    
    fl.rst = open(suffix * "/" * "irt.rst" * name_end, "w+")
    fl.res = open(suffix * "/" * "irt.res" * name_end, "w")
    fl.col = open(suffix * "/" * "irt.col" * name_end, "w")
    fl.sex = open(suffix * "/" * "irt.sex" * name_end, "w")
    fl.trj = open(suffix * "/" * "irt.trj" * name_end, "w")
    fl.wsk = open(suffix * "/" * "irt.wsk" * name_end, "w")
    fl.tpp = open(suffix * "/" * "irt.tpp" * name_end, "w")
    fl.bur = open(suffix * "/" * "irt.bur" * name_end, "w")

end #open_files

# This function closes files.
function close_files(fl)

    close(fl.rst)
    close(fl.res)
    close(fl.col)
    close(fl.sex)
    close(fl.trj)
    close(fl.wsk)
    close(fl.tpp)
    close(fl.bur)

end #close_files

# This function writes average values as functions of the parameters.
function write_avr(sval, svala, av, icall, fl)

    println(fl.out, "icall=", icall)
    ivar = 0
    icall == 1 ? println(fl.out, "\navr var") : nothing
        
    if sval.scan_type == 't'
        printfmt(fl.avr, "{1:f} {2:d} ", svala.par, svala.ipar)
	if icall == 1
	    ivar += 1
	    printfmtln(fl.out, "{1:d} svala.par", ivar)
	    printfmtln(fl.avc, "{1:d} para Float64[]", ivar)
	    ivar += 1
	    printfmtln(fl.out, "{1:d} svala.ipar", ivar)
	    printfmtln(fl.avc, "{1:d} ipara Int64[]", ivar)
        end
    end

    if sval.scan_type in Set(['y','t'])
        printfmt(fl.avr, "{1:f} {2:d} {3:d}", sval.par, sval.ipar, sval.irepeat)
 	if icall == 1
	    ivar += 1
	    printfmtln(fl.out, "{1:d} sval.par", ivar)
	    printfmtln(fl.avc, "{1:d} par Float64[]", ivar)
	    ivar += 1
	    printfmtln(fl.out, "{1:d} sval.ipar", ivar)
	    printfmtln(fl.avc, "{1:d} ipar Int64[]", ivar)
	    ivar += 1
	    printfmtln(fl.out, "{1:d} sval.irepeat", ivar)
	    printfmtln(fl.avc, "{1:d} irepeat Int64[]", ivar)
        end
   end
      
    for ipop in 1:av[1].Npop
        printfmt(fl.avr, " {1:f} {2:f} {3:f}", av[ipop].V_avt,
	  av[ipop].V_avt_spk_clamp, av[ipop].chi)
 	if icall == 1
	    ivar += 1
	    printfmtln(fl.out, "{1:d} V_avt[{2:d}]", ivar, ipop)
	    printfmtln(fl.avc, "{1:d} V_avt{2:d} Float64[]", ivar, ipop)
	    ivar += 1
	    printfmtln(fl.out, "{1:d} V_avt_spk_clamp[{2:d}]", ivar, ipop)
	    printfmtln(fl.avc, "{1:d} V_avt_spk_clamp{2:d} Float64[]", ivar,
	      ipop)
	    ivar += 1
	    printfmtln(fl.out, "{1:d} chi[{2:d}]", ivar, ipop)
	    printfmtln(fl.avc, "{1:d} chi{2:d} Float64[]", ivar, ipop)
	end
	
	printfmt(fl.avr, " {1:d} {2:f} {3:f} {4:f}", av[ipop].n_silent,
	  av[ipop].CV_pop_av, av[ipop].CVt_pop_av, av[ipop].T_pop_av)
 	if icall == 1
	    ivar += 1
	    printfmtln(fl.out, "{1:d} n_silent[{2:d}]", ivar, ipop)
	    printfmtln(fl.avc, "{1:d} n_silent{2:d} Float64[]", ivar, ipop)
	    ivar += 1
	    printfmtln(fl.out, "{1:d} CV_pop_av[{2:d}]", ivar, ipop)
	    printfmtln(fl.avc, "{1:d} CV_pop_av{2:d} Float64[]", ivar, ipop)
	    ivar += 1
	    printfmtln(fl.out, "{1:d} CVt_pop_av[{2:d}]", ivar, ipop)
	    printfmtln(fl.avc, "{1:d} CVt_pop_av{2:d} Float64[]", ivar, ipop)
	    ivar += 1
	    printfmtln(fl.out, "{1:d} T_pop_av[{2:d}]", ivar, ipop)
	    printfmtln(fl.avc, "{1:d} T_pop_av{2:d} Float64[]", ivar, ipop)
	end
	
	printfmt(fl.avr, " {1:f}", av[ipop].Mfr)
 	if icall == 1
	    ivar += 1
	    printfmtln(fl.out, "{1:d} Mfr[{2:d}]", ivar, ipop)
	    printfmtln(fl.avc, "{1:d} Mfr{2:d} Float64[]", ivar, ipop)
	end

        if ipop <= 2
	    printfmt(fl.avr, " {1:f} {2:f}", av[ipop].CV_bur_pop_av,
	      av[ipop].CVt_bur_pop_av)
 	    if icall == 1
	        ivar += 1
	        printfmtln(fl.out, "{1:d} CV_bur_pop_av[{2:d}]", ivar, ipop)
	        printfmtln(fl.avc, "{1:d} CV_bur_pop_av{2:d} Float64[]", ivar,
		  ipop)
	        ivar += 1
	        printfmtln(fl.out, "{1:d} CVt_bur_pop_av[{2:d}]", ivar, ipop)
	        printfmtln(fl.avc, "{1:d} CVt_bur_pop_av{2:d} Float64[]", ivar,
		  ipop)
	    end
	end
    end

#   0 instead of av[1].mode
##    printfmt(fl.avr, " {1:d} {2:f} {3:f}", 0, av[1].Tper,
##      av[1].Tper_from_cc_av)
      
    printfmt(fl.avr, " {1:d} {2:f} {3:f}", av[1].mode, av[1].Tper,
      av[1].Tper_from_cc_av)
      
    if icall == 1
	ivar += 1
	printfmtln(fl.out, "{1:d} mode", ivar)
	printfmtln(fl.avc, "{1:d} mode Int64[]", ivar)
	ivar += 1
	printfmtln(fl.out, "{1:d} Tper", ivar)
	printfmtln(fl.avc, "{1:d} Tper Float64[]", ivar)
	ivar += 1
	printfmtln(fl.out, "{1:d} Tper_from_cc", ivar)
	printfmtln(fl.avc, "{1:d} Tper_from_cc Float64[]", ivar)
    end

    printfmt(fl.avr, " {1:f} {2:f} {3:f} {4:f}", av[1].Sp_av, av[1].Sp_sd, av[1].Sm_av, av[1].Sm_sd)
    
    if icall == 1
	ivar += 1
	printfmtln(fl.out, "{1:d} Sp_av", ivar)
	printfmtln(fl.avc, "{1:d} Sp_av Float64[]", ivar)
	ivar += 1
	printfmtln(fl.out, "{1:d} Sp_sd", ivar)
	printfmtln(fl.avc, "{1:d} Sp_sd Float64[]", ivar)
	ivar += 1
	printfmtln(fl.out, "{1:d} Sm_av", ivar)
	printfmtln(fl.avc, "{1:d} Sm_av Float64[]", ivar)
	ivar += 1
	printfmtln(fl.out, "{1:d} Sm_sd", ivar)
	printfmtln(fl.avc, "{1:d} Sm_sd Float64[]", ivar)
    end

    for ipop in 1:av[1].Npop
        printfmt(fl.avr, " {1:f}", av[ipop].spk_during_prebot_avt_avpop)
	
  	if icall == 1
	    ivar += 1
	    printfmtln(fl.out, "{1:d} spk_du[{2:d}]", ivar, ipop)
	    printfmtln(fl.avc, "{1:d} spk_du{2:d} Float64[]", ivar, ipop)
	end
   end

    printfmt(fl.out, "{1:c}", '\n')
    printfmt(fl.avr, "{1:c}", '\n')
    flush(fl.avr)   
    flush(fl.avc)   

    sval.scan_type == 't' ? printfmt(fl.aws, "{1:f} {2:d} ",
      svala.par, svala.ipar) : nothing
    sval.scan_type in Set(['y','t']) ? printfmt(fl.aws, "{1:f} {2:d} {3:d} ",
      sval.par, sval.ipar, sval.irepeat) : nothing
end #write_avr

# This function writes peak heights of consecuite whisking cycles as
# functions of the parameters.
function write_hei(sval, svala, height_maxima, del_time_maxima, fl)

    sval.scan_type == 't' ? printfmt(fl.hei, "{1:f} {2:d} ",
      svala.par, svala.ipar) : nothing
    sval.scan_type in Set(['y','t']) ? printfmt(fl.hei, "{1:f} {2:d} {3:d}",
      sval.par, sval.ipar, sval.irepeat) : nothing
    for peak in height_maxima
        printfmt(fl.hei, " {1:f}", peak)
    end
    printfmt(fl.hei, "{1:c}", '\n')
    
    flush(fl.hei)

    sval.scan_type == 't' ? printfmt(fl.dtm, "{1:f} {2:d} ",
      svala.par, svala.ipar) : nothing
    sval.scan_type in Set(['y','t']) ? printfmt(fl.dtm, "{1:f} {2:d} {3:d}",
      sval.par, sval.ipar, sval.irepeat) : nothing
    for det_time in del_time_maxima
        printfmt(fl.dtm, " {1:f}", det_time)
    end
    printfmt(fl.dtm, "{1:c}", '\n')
    
    flush(fl.dtm)

end #write_hei

# This function call the fucntions for computing hei if mode == 2.	    
function call_hei(rstpbstr, sval, partohop, suffix, runpar, av, fl)

    height_maxima, del_time_maxima = find_whisking_amplitudes(rstpbstr, sval,
      partohop, suffix, av[1].mode, fl)
    write_hei(sval, svala, height_maxima, del_time_maxima, fl)
    println("height_maxima=", height_maxima)
    println(" del_time_maxima=", del_time_maxima)
    
end #call_hei

# main

dir_dat = pwd() * "/" # "../dat5/"
include("irt_types.jl")
include("irt_sim.jl")
include("hop_cal.jl")
include("hop_types.jl")
include("bur_cal.jl")

if length(ARGS) >= 1
    suffix = ARGS[1]
else
    suffix = "a1"
end

fl = fl_st()
!isdir(suffix) ? mkdir(suffix) : nothing
fl.avr = open(suffix * "/" * "irt.avr", "w")
fl.avc = open(suffix * "/" * "irt.avc", "w")
fl.hei = open(suffix * "/" * "irt.hei", "w")
fl.dtm = open(suffix * "/" * "irt.dtm", "w")
fl.aws = open(suffix * "/" * "irt.aws", "w")
fl.out = open(suffix * "/" * "irt.out", "w")

println("suffix=", suffix)
include(dir_dat * suffix * "/irt_const.jl")

Random.seed!(sval.seed)

av = Array{avr_val,1}(undef, 0)  #(undef, 5)
rstpbstr = rst_pb_str()
partohop = par_to_hop()

if sval.scan_type == 'n'
    runpar.open_file_for_par = true
    open_files(sval.scan_type, suffix, 0, 0, 1, fl)
#   rst_ar = constract_rst_ar(netpar.Npop)
    one_par(suffix, netpar, sval, av, partohop, rstpbstr, fl)
    write_avr(sval, svala, av, 1, fl)
    runpar.whisk_cal ? call_hei(rstpbstr, sval, partohop,suffix, runpar, av,
      fl) : av[1].mode = 0
    close_files(fl)
elseif sval.scan_type == 'y'
    runpar.open_file_for_par = sval.open_file_for_par
    for sval.ipar in 0:sval.npar
        sval.par = sval.npar == 0 ? sval.parmin :  sval.parmin +
	  (sval.parmax - sval.parmin) * sval.ipar / sval.npar
	println("par=", sval.par, " ipar=", sval.ipar)
	modify_parameter(sval)
	for sval.irepeat in 1:sval.nrepeat
	    println("irepeat=", sval.irepeat)
	    sval.open_file_for_par ? open_files(sval.scan_type, suffix,
	      0, sval.ipar, sval.irepeat, fl) : nothing
#           rst_ar = constract_rst_ar(netpar.Npop)
            one_par(suffix, netpar, sval, av, partohop, rstpbstr, fl)
	    write_avr(sval, svala, av, sval.ipar * sval.nrepeat + sval.irepeat,
	      fl)
            runpar.whisk_cal ? call_hei(rstpbstr, sval, partohop,suffix,
	      runpar, av, fl) :  av[1].mode = 0
	    sval.open_file_for_par ? close_files(fl) : nothing
	end
    end
elseif sval.scan_type == 't'
    runpar.open_file_for_par = sval.open_file_for_par
    for svala.ipar in 0:svala.npar, sval.ipar in 0:sval.npar
        svala.par = svala.npar == 0 ? svala.parmin :  svala.parmin +
	  (svala.parmax - svala.parmin) * svala.ipar / svala.npar
        sval.par = sval.npar == 0 ? sval.parmin :  sval.parmin +
	  (sval.parmax - sval.parmin) * sval.ipar / sval.npar
	println("para=", svala.par, " ipara=", svala.ipar, " par=", sval.par,
	  " ipar=", sval.ipar)
	modify_parameter(svala)
	modify_parameter(sval)
	for sval.irepeat in 1:sval.nrepeat
	    println("irepeat=", sval.irepeat)
	    sval.open_file_for_par ? open_files(sval.scan_type, suffix,
	      svala.ipar, sval.ipar, sval.irepeat, fl) : nothing
#           rst_ar = constract_rst_ar(netpar.Npop)
            one_par(suffix, netpar, sval, av, partohop, rstpbstr, fl)
	    write_avr(sval, svala, av, svala.ipar * sval.npar * sval.nrepeat +
	      sval.ipar * sval.nrepeat + sval.irepeat, fl)
            runpar.whisk_cal ? call_hei(rstpbstr, sval, partohop,suffix, runpar,
	      av, fl) :  av[1].mode = 0
	    sval.open_file_for_par ? close_files(fl) : nothing
	end
    end
else
    println("wrong scan_type=", sval.scan_type)
    exit(0)
end

close(fl.avr)
close(fl.avc)
close(fl.hei)
close(fl.dtm)
close(fl.aws)
close(fl.out)
