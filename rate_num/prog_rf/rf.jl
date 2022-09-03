# This program simulates the differential equations of the rate model.
# preBot->vIRtr->FN

using Formatting
using Printf
using DataStructures

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

# This function opens file for writing.
function open_files(suffix, strfl, fl)
    fl.col = open(suffix * "/rf.col" * strfl, "w")
    fl.grb = open(suffix * "/rf.grb" * strfl, "w")
    fl.pek = open(suffix * "/rf.pek" * strfl, "w")
end #open_files

# This program closes files for writing.
function close_files(fl)
    close(fl.col)
    close(fl.grb)
    close(fl.pek)
end #close_files

# This function writes average values as functions of the parameters.
function write_avr(sval, svala, av, fl)

    sval.scan_type == 't' ? printfmt(fl.avr, "{1:f} {2:d} ",
      svala.par, svala.ipar) : nothing
    sval.scan_type in Set(['y','t']) ? printfmt(fl.avr, "{1:f} {2:d} {3:d}",
      sval.par, sval.ipar, sval.irepeat) : nothing
    [printfmt(fl.avr, " {1:f} {2:f}", av.Mav[ipop], av.Msd[ipop]) for
      ipop in 1:length(av.Mav)]
    printfmt(fl.avr, " {1:f}", av.Tper)
    printfmt(fl.avr, " {1:f} {2:f}", av.Mup[1], av.Mup[3])
    printfmt(fl.avr, " {1:f} {2:f} {3:f} {4:f}", av.theta_av, av.theta_sd,
      av.whisk_amp_avr, av.whisk_setup)
    printfmt(fl.avr, "\n")
    flush(fl.avr)
   
end #write_avr

# main

dir_dat = pwd() * "/" #"../dat/"
include("rf_types.jl")
include("rf_sim.jl")
include("rf_peaks.jl")

if length(ARGS) >= 1
    suffix = ARGS[1]
else
    suffix = "a1"
end

fl = fl_st()
!isdir(suffix) ? mkdir(suffix) : nothing
fl.avr = open(suffix * "/" * "rf.avr", "w")
fl.wsk = open(suffix * "/" * "rf.wsk", "w")
fl.tim = open(suffix * "/" * "rf.tim", "w")
fl.out = open(suffix * "/" * "rf.out", "w")

include(dir_dat * suffix * "/rf_const.jl")

av = avr_val()

if sval.scan_type == 'n'
    open_files(suffix, "", fl)
    one_par(suffix, netpar, sval, av, fl)
    write_avr(sval, svala, av, fl)
    close_files(fl)
elseif sval.scan_type == 'y'
    for sval.ipar in 0:sval.npar
        sval.par = sval.npar == 0 ? sval.parmin :  sval.parmin +
	  (sval.parmax - sval.parmin) * sval.ipar / sval.npar
	println("par=", sval.par, " ipar=", sval.ipar)
      	modify_parameter(sval)
	for sval.irepeat in 1:sval.nrepeat
	    println("irepeat=", sval.irepeat)
	    runpar.open_file_for_par ? open_files(suffix,
	      ".$(sval.ipar).$(sval.irepeat)", fl) : nothing
	    printfmtln(fl.wsk, "{1:f}", sval.par)
	    printfmtln(fl.tim, "{1:f}", sval.par)
            one_par(suffix, netpar, sval, av, fl)
            write_avr(sval, svala, av, fl)
            runpar.open_file_for_par ? close_files(fl) : nothing
        end
    end
elseif sval.scan_type == 't'
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
            runpar.open_file_for_par ? open_files(suffix,
	      ".$(svala.ipar).$(sval.ipar).$(sval.irepeat)", fl) : nothing
	    printfmtln(fl.wsk, "{1:f} {2:f}", svala.par, sval.par)
 	    printfmtln(fl.tim, "{1:f} {2:f}", svala.par, sval.par)
           one_par(suffix, netpar, sval, av, fl)
            write_avr(sval, svala, av, fl)
            runpar.open_file_for_par ? close_files(fl) : nothing
        end
    end
else
    println("wrong scan_type=", sval.scan_type)
    exit(0)
end

close(fl.avr)
close(fl.wsk)
close(fl.tim)
close(fl.out)
