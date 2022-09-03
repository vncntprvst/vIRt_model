# Two-vIRT system. Fast-slow analysis.
# Computing the oscillation time period.

using Formatting
using Printf
using DataStructures
#using Statistics
#using DSP
using NLsolve

#mutable struct net_par
#    I0::Float64
#    taua::Float64
#    taus::Float64
#    Ja::Float64
#    Jinter::Float64
#    Jintra::Float64
#    beta::Float64
#    gamma::Float64
#    gKm::Float64
#    DelV::Float64
#    gintra::Float64
#    ginter::Float64
        
#    net_par() = new()
#end

#mutable struct run_par
#    ncut::Int64
#    eps_from_bif::Float64
#    parmax::Float64
        
#    run_par() = new()
#end

mutable struct fl_st_tvnc
    res::IOStream
    mfr::IOStream
    amp::IOStream
    set::IOStream
    out::IOStream
    
    fl_st_tvnc() = new()
end

# This function substitutes the parameter value xval in the appropriate
# parameter value.
function substitute_parameter(xval, netpar, fl)

	if sval.parname == ["netpar.ginter"]
            netpar.Jinter = xval
	elseif sval.parname == ["netpar.I0"]
	    netpar.I0 = xval
	else
	    println("parname=", sval.parname, " is not in the list!")
	    exit(0)
	end

end #substitute_parameter

# This function finds the bifurcation points
function compute_bif(netpar, fl)

    JpTr = (1 / netpar.beta) * (1/netpar.taus + 1/netpar.taua +
      netpar.beta * netpar.Ja / netpar.taua)
    JpDet = (1 / netpar.beta) * (1/netpar.taus) * (1 + netpar.beta * netpar.Ja)

    return [netpar.Jintra + JpTr, netpar.Jintra + JpDet]
end #compute_bif

# This function computes M2cal for a specific value of TT.
function compute_m2cal!(m2calar, TTar, netpar, fl)

    TT = TTar[1]
    
    Jtilde = netpar.Ja * netpar.beta / (1.0 + netpar.taus * netpar.beta * netpar.Jintra)

    fracIJ = (netpar.I0 - netpar.I0_threshold) * Jtilde / (1.0 + Jtilde)
    expJtT = exp(-(1 + Jtilde) * TT / netpar.taua)
    expTt  = exp(-TT / netpar.taua)
    a0 = fracIJ * (1 - expJtT) * expTt / (1 - expJtT * expTt)
    aT = a0 * expJtT + fracIJ * (1.0 - expJtT)

    s1 = netpar.taus * netpar.beta * (netpar.I0 - netpar.I0_threshold - aT) /
      (1 + netpar.taus * netpar.beta * netpar.Jintra) 
    m2cal = netpar.I0 - netpar.I0_threshold - netpar.Jinter * s1 - a0

    m2calar[1] = m2cal

end #compute_m2cal!

# This function computes the half-time-period TT.
function compute_TT(TTar_init, netpar, fl)

#   TTar = [200.0]
    m2calar = [0.0]

    compute_m2cal_modified!(m2calar, TTar) =
      compute_m2cal!(m2calar, TTar, netpar, fl)
    sol_m2cal = nlsolve(compute_m2cal_modified!, TTar_init, autodiff = :forward)
    TTar = sol_m2cal.zero
#   println("TTar=", TTar, " sol_m2cal=", sol_m2cal)
    
    return TTar[1]
end #compute_TT

# This function find TT for a parameter interval.
function find_TT_for_par(bif_values, netpar, fl)

    netpar.Jinter = (netpar.DelV / netpar.taus) * netpar.ginter

    if sval.parname == ["netpar.ginter"]
        xval_min = bif_values[2] - runpar.eps_from_bif
        xval_max = bif_values[1] + runpar.eps_from_bif
    elseif sval.parname == ["netpar.I0"]
        xval_min = netpar.I0_min
        xval_max = netpar.I0_max
    end

    if sval.parname == ["netpar.I0"] &&
      (netpar.Jinter < bif_values[1] + runpar.eps_from_bif ||
       netpar.Jinter > bif_values[2] - runpar.eps_from_bif)
        println("Jinter=", netpar.Jinter, " is outside of the range: [",
	  xval_min, ",", xval_max, "]")
	return
    end

    TTar_init = [100.0]
    for icut in 0:runpar.ncut
        xval = xval_min + icut * (xval_max - xval_min) / runpar.ncut
        substitute_parameter(xval, netpar, fl)

        if netpar.I0 > netpar.I0_threshold
            TT = compute_TT(TTar_init, netpar, fl)

            if sval.parname == ["netpar.ginter"]
                printfmtln(fl.res, "{1:f} {2:f} {3:f}", netpar.Jinter *
	          (netpar.taus / netpar.DelV), TT, netpar.Jinter)
	    elseif sval.parname == ["netpar.I0"]
	        printfmtln(fl.res, "{1:f} {2:f}", netpar.I0, TT)
	    end
	    TTar_init[1] = TT
	end
    end
    
end #find_TT_for_par(netpar, fl)

# This function computes analytically the integral over Mfr for the
# periodic state.
function compute_Mfr_integral(TT, netpar, fl)

    Jtilde = netpar.Ja * netpar.beta / (1.0 + netpar.taus * netpar.beta *
      netpar.Jintra)
    fracIJ = (netpar.I0 - netpar.I0_threshold) * Jtilde / (1.0 + Jtilde)
    expJtT = exp(-(1 + Jtilde) * TT / netpar.taua)
    expTt  = exp(-TT / netpar.taua)
	
    factor = netpar.beta / (1.0 + netpar.taus * netpar.beta * netpar.Jintra)
	
	Int1 = 0.5 * factor * (netpar.I0 - netpar.I0_threshold - fracIJ)
        Expf1 = exp(TT / netpar.taua) - 1.0
	Expf2 = exp(TT / netpar.taua) - exp(-(1.0 + Jtilde) * TT / netpar.taua)
	Expf3 = 1.0 - exp(-(1.0 + Jtilde) * TT / netpar.taua)
	Int21 = (0.5 * factor / TT) * fracIJ * (netpar.taua / (1.0 + Jtilde)) 
	Int2 = Int21 * (Expf1 / Expf2) * Expf3
	Mfr_int = Int1 + Int2
#	printfmtln("TT={1:f}", TT)
#	printfmtln("Int1={1:f} Int21={2:f} Expf1={3:f} Expf2={4:f} " *
#         "Expf3={5:f} Int2={6:f} Mfr_int={7:f}", Int1, Int21, Expf1, Expf2,
#	  Expf3, Int2, Mfr_int)
	
# Check manuscript
    Tv = 2.0 * TT
    Itilde = netpar.I0 - netpar.I0_threshold

    A1 = netpar.beta * Itilde /
      (2.0 * (1.0 + netpar.taus * netpar.beta * netpar.Jintra) * (1 + Jtilde))
    Be1 = exp(Tv / (2 * netpar.taua)) - 1.0
    Be2 = 1.0 - exp(-(1 + Jtilde) * Tv / (2 * netpar.taua))
    Be3 = exp(Tv / (2 * netpar.taua)) -
      exp(-(1 + Jtilde) * Tv / (2 * netpar.taua))
    A21 = (netpar.beta * Itilde * Jtilde  * netpar.taua /
      (Tv * (1.0 + netpar.taus * netpar.beta * netpar.Jintra) *
      (1 + Jtilde)^2)) 
    A2 = A21 * (Be1 * Be2 / Be3)
    M_check = A1 + A2
#   printfmtln("A1={1:f} A21={2:f} Be1={3:f} Be2={4:f} Be3={5:f} A2={6:f} " *
#     "M_check={7:f}", A1, A21, Be1, Be2, Be3, A2, M_check)

# Check another version of the manuscript
    expTt = exp(-TT / netpar.taua)
    expJtT = exp(-(1 + Jtilde) * TT / netpar.taua)
    expJtTT = exp(-(2 + Jtilde) * TT / netpar.taua)
    
    factor_c1 = netpar.beta * Itilde /
      ((1.0 + netpar.taus * netpar.beta * netpar.Jintra) * (1 + Jtilde))

    numer = Jtilde * netpar.taua * (1.0 - expTt) * (1.0  - expJtT)
    denom = Tv * (1.0 + Jtilde) * (1.0 - expJtTT)

    M_check_A = factor_c1 * (0.5 + (numer / denom))


    factor_c1_B = Jtilde * Itilde / (2 * (1 + Jtilde) * netpar.Ja)
    numer_B = 2 * Jtilde * netpar.taua * (1.0 - expTt) * (1.0  - expJtT)
    denom_B = (1.0 + Jtilde) * Tv * (1.0 - expJtTT)
    
    M_check_B = factor_c1_B * (1 + (numer_B / denom_B))

#   printfmtln("Check: TT={1:f} Mfr_int={2:f} M_check={3:f} M_check_A={4:f} " *
#     "M_check_B={4:f}", TT, Mfr_int, M_check, M_check_A)
 
    return Mfr_int
end #compute_Mfr_integral

# This function computes the firing rates of the sub-nuclei.
function compute_Mfr_for_par(bif_values, netpar, fl)

    ncut_m = runpar.ncut / 10
    println("bif_values=", bif_values)
    
#   Two populations fire tonically.
    if sval.parname == ["netpar.ginter"]
        for icut in 0:ncut_m
            xval = bif_values[1] * icut / ncut_m
            substitute_parameter(xval, netpar, fl)
	    ginter = (netpar.taus / netpar.DelV) * netpar.Jinter
	
	    Mfr = 1000.0 * netpar.beta * (netpar.I0 - netpar.I0_threshold) /
	      (1.0 + netpar.beta * (netpar.taus * (netpar.Jintra +
	      netpar.Jinter) + netpar.Ja))
	    printfmtln(fl.mfr, "{1:f} {2:f}", ginter, Mfr)
#	    printfmtln("{1:f} {2:f} {3:f} {4:f} {5:f}", ginter, Mfr,
#	      netpar.Jinter, netpar.Jintra, netpar.Ja)
	end
    end

#   Limit cycle. The two populations fire alternately.
    TTar_init = [30.0]
    Mfr_last = 0.0
    
    if sval.parname == ["netpar.ginter"]
       xval_min = bif_values[1]
       xval_max = bif_values[2]
    elseif sval.parname == ["netpar.I0"]
       xval_min = netpar.I0_min
       xval_max = netpar.I0_max
    end
    
    for icut in 0:ncut_m
        xval = xval_min + (xval_max - xval_min) * icut / ncut_m
        substitute_parameter(xval, netpar, fl)
	ginter = (netpar.taus / netpar.DelV) * netpar.Jinter

        if netpar.I0 > netpar.I0_threshold
            TT = compute_TT(TTar_init, netpar, fl)
    	    TTar_init[1] = TT
	    Mfr = compute_Mfr_integral(TT, netpar, fl)
	    x_write = sval.parname == ["netpar.ginter"] ? ginter : xval
	    printfmtln(fl.mfr, "{1:f} {2:f} {3:f}", x_write, 1000.0 * Mfr, TT)
	    Mfr_last = Mfr
	end
    end

#   One population fires tonically, one is silent.    
    if sval.parname == ["netpar.ginter"]
        ginter_bif = bif_values[2] * (netpar.taus / netpar.DelV)

        Jinter_max = runpar.parmax * (netpar.DelV / netpar.taus)
        for icut in 0:ncut_m
            xval = bif_values[2] + (Jinter_max - bif_values[2]) * icut / ncut_m
            substitute_parameter(xval, netpar, fl)
	    ginter = (netpar.taus / netpar.DelV) * netpar.Jinter
	    
	    if netpar.I0 > netpar.I0_threshold
	        Mfr = 1000.0 * netpar.beta * (netpar.I0 - netpar.I0_threshold) /
	          (1.0 + netpar.beta * (netpar.taus * netpar.Jintra +
		  netpar.Ja))
	        printfmtln(fl.mfr, "{1:f} {2:f}", ginter, Mfr)
            end
	end
	
        printfmtln(fl.mfr, "  ")
        printfmtln(fl.mfr, "{1:f} {2:f}", ginter_bif, 1000.0 * Mfr_last)
        printfmtln(fl.mfr, "{1:f} {2:f}", ginter_bif, 0.0)
        printfmtln(fl.mfr, "{1:f} {2:f}", runpar.parmax, 0.0)
    end

end #compute_Mfr_for_par

# This function computes analytically vlaues for the whisking amplitude and
# set point.
function compute_wsk_amp_set_for_par(bif_values, netpar, fl)

#   Amplitude
    println(fl.amp, "-100.0 0.0")
    printfmtln(fl.amp, "{1:f} 0.0", bif_values[1] *
      (netpar.taus / netpar.DelV) - netpar.gintra)
    println(fl.amp, "  ")
    printfmtln(fl.amp, "{1:f} 0.0", bif_values[2] *
      (netpar.taus / netpar.DelV) - netpar.gintra)
    println(fl.amp, "100.0 0.0")

    ncut_m = runpar.ncut / 10
    
#   Set point

#   Two populations fire tonically.
    for icut in 0:ncut_m
        xval = bif_values[1] * icut / ncut_m
        substitute_parameter(xval, netpar, fl)
        ginter = (netpar.taus / netpar.DelV) * netpar.Jinter
    
        Mr = netpar.beta * (netpar.I0 - netpar.I0_threshold) /
          (1.0 + netpar.beta * (netpar.taus * (netpar.Jintra +
          netpar.Jinter) + netpar.Ja))

        MF = (netpar.beta_F / (1.0 + netpar.Ja_F * netpar.beta_F)) *
	  (netpar.I0_F - netpar.I0_threshold_F - netpar.JFr * netpar.taus * Mr)
	MF < 0.0 ? MF = 0.0 : nothing

        force_term = force_from_MF(MF, netpar, runpar, fl)

        theta_set = netpar.tauw * force_term 

#       printfmtln("{1:f} {2:f} {3:f}", ginter - netpar.gintra, theta_set, Mr,
#       MF, force_term)
        printfmtln(fl.set, "{1:f} {2:f} {3:f} {4:f} {5:f}",
	  ginter - netpar.gintra, theta_set, Mr, MF, force_term)
    end

    println(fl.set, "  ")

#   One population fires tonically, one is silent.    
    Jinter_max = 100.0 * (netpar.DelV / netpar.taus)
    for icut in 0:ncut_m
        xval = bif_values[2] + (Jinter_max - bif_values[2]) * icut / ncut_m
        substitute_parameter(xval, netpar, fl)
        ginter = (netpar.taus / netpar.DelV) * netpar.Jinter
        
        if netpar.I0 > netpar.I0_threshold
            Mr = netpar.beta * (netpar.I0 - netpar.I0_threshold) /
              (1.0 + netpar.beta * (netpar.taus * netpar.Jintra + netpar.Ja))
	      
            MF = (netpar.beta_F / (1.0 + netpar.Ja_F * netpar.beta_F)) *
	      (netpar.I0_F - netpar.I0_threshold_F)
	    MF < 0.0 ? MF = 0.0 : nothing

            force_term = force_from_MF(MF, netpar, runpar, fl)

            theta_set = netpar.tauw * force_term 

            printfmtln(fl.set, "{1:f} {2:f} {3:f} {4:f} {5:f}",
	      ginter - netpar.gintra, theta_set, Mr, MF, force_term)
        end
    end
    
    println(fl.set, "  ")
    printfmtln(fl.set, "{1:f} 0.0", bif_values[2] *
      (netpar.taus / netpar.DelV) - netpar.gintra)
    println(fl.set, "100.0 0.0")

end #compute_wsk_amp_set_for_par

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

# main

if length(ARGS) >= 1
    suffix = ARGS[1]
else
    suffix = "no_pb_x_ginter_intra" #"a5" "a4" #"a3"
end

current_dir = pwd()

println("suffix=", suffix)
fl = fl_st_tvnc()
fl.res = open(suffix * "/tvnc.res", "w")
fl.mfr = open(suffix * "/tvnc.mfr", "w")
fl.amp = open(suffix * "/tvnc.amp", "w")
fl.set = open(suffix * "/tvnc.set", "w")
fl.out = open(suffix * "/tvnc.out", "w")

# reference parameter values
#include("../rate_num/prog/numr_types.jl")
#include("../prog_rf/rf_types.jl")
include("../rate_num/prog_rf/rf_types.jl")
netpar = net_par()
runpar = run_par()
#include("../rate_num/datfig/" * suffix * "/numr_const.jl")
include(current_dir * "/" * suffix * "/" * "/rf_const.jl")

netpar.I0_min = 0.0
netpar.I0_max = 40.0

netpar.Ja = 1000.0 * netpar.gamma * netpar.gKm
netpar.Jintra = netpar.gintra * netpar.DelV / netpar.taus
netpar.Jinter = netpar.ginter * netpar.DelV / netpar.taus
printfmtln(fl.out, "Ja={1:f} Jintra={2:f} Jinter={3:f}", netpar.Ja,
  netpar.Jintra, netpar.Jinter) 

netpar.Ja_F = 1000.0 * netpar.gamma_F * netpar.gKm_F 
netpar.JFr = netpar.gFr * netpar.DelV_F / netpar.taus
printfmtln(fl.out, "Ja_F={1:f} JFr={2:f}", netpar.Ja_F, netpar.JFr)


runpar.ncut = 1000
runpar.eps_from_bif = 0.001
runpar.parmax = 55.0
bif_values = compute_bif(netpar, fl)
printfmtln("bif_values={1:f} {2:f}", bif_values[1], bif_values[2])

printfmtln(fl.out, "DelV={1:f} bif_values={2:f} {3:f}", netpar.DelV,
  bif_values[1] * (netpar.taus / netpar.DelV), bif_values[2] *
  (netpar.taus / netpar.DelV))

find_TT_for_par(bif_values, netpar, fl)

compute_Mfr_for_par(bif_values, netpar, fl)

if sval.parname == ["netpar.ginter"]
    compute_wsk_amp_set_for_par(bif_values, netpar, fl)
else
    println(fl.amp, "  ")
end

close(fl.res)
close(fl.mfr)
close(fl.amp)
close(fl.set)
close(fl.out)
