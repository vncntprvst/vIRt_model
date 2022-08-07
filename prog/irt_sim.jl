# This file contains the functions for simulating the dynamical system.

#include("irt_par.jl")

# This function defines the sigmoid function Gammaf
function Gammaf(VV, theta, sigma)
    return 1.0/(1.0+exp(-(VV-theta)/sigma))
end #Gammaf

#This function calculates the steady-state variables for a specific V
#for IRT neurons.
function steady_state_var_I_br(Varc_one, ion, cellpar, fl)
    Vc = Varc_one[1]

    ninf  = 1 / (1+exp((Vc-cellpar.thetan) /cellpar.sigman))
    ksinf = 1 / (1+exp((Vc-cellpar.thetaks)/cellpar.sigmaks))
    
    FF = 1 / (1 + exp(-(Vc-cellpar.theta_syn)/cellpar.sigma_syn))

    Varc_one[2] = ninf
    Varc_one[3] = ksinf
    Varc_one[4] = (cellpar.alpha * FF) / (cellpar.alpha * FF + cellpar.beta)
    
end #steady_state_var_I_br

#This function updates the variables of an IRT neuron model.
function update_cell_I_br!(Varc_one, kout_one, ion, cellpar, Isyn_cont,
         Iel_cont, s_ext, netpar, it, time, fl)

    Vc = Varc_one[1]
    nc = Varc_one[2]
    ksc = Varc_one[3]
    sync = Varc_one[4]

    minf  = 1 / (1+exp((Vc-cellpar.thetam) /cellpar.sigmam)) 
    ninf  = 1 / (1+exp((Vc-cellpar.thetan) /cellpar.sigman))
    pinf  = 1 / (1+exp((Vc-cellpar.thetap) /cellpar.sigmap))
    ksinf = 1 / (1+exp((Vc-cellpar.thetaks)/cellpar.sigmaks))

    taun  = cellpar.taunb0 +
      cellpar.taunb/cosh((Vc-cellpar.thetan)/(2*cellpar.sigman))
    tauks = cellpar.tauksb0 +
      cellpar.tauksb/cosh((Vc-cellpar.thetaks)/(2*cellpar.sigmaks))

    INa  = cellpar.gNa  * minf*minf*minf * (1-nc) * (Vc-cellpar.VNa)
    INap = cellpar.gNap * pinf * (Vc-cellpar.VNa)
    IKdr = cellpar.gKdr[ion] * nc*nc*nc*nc * (Vc-cellpar.VK)
    IKs  = cellpar.gKs  * ksc * (Vc-cellpar.VK)
    IL   = netpar.gLar[ion] * (Vc-cellpar.VL)
    
    FF = 1 / (1 + exp(-(Vc-cellpar.theta_syn)/cellpar.sigma_syn))

    Isyn_ext = netpar.gsyn_ext_ar[ion] * s_ext * (Vc - cellpar.Vsyn_ext)
#   ion == 1 ? println(fl.col, time, " Ise=",  Isyn_ext) : nothing

#V
    kout_one[1] = (-INa - INap - IKdr - IKs - IL - Isyn_cont - Iel_cont -
      Isyn_ext + cellpar.Iapp) / cellpar.Cms
#    ion == 1 ? printfmtln(fl.col,
#      "{1:d} {2:f} {3:f} {4:f} {5:f} {6:f} {7:f} {8:f}", it,  kout_one[1],
#      INa, INap, IKdr, IKs, IL, cellpar.Iapp) : nothing
#n
    kout_one[2] = cellpar.phi*(ninf-nc)/taun
#ks
    kout_one[3] = (ksinf-ksc)/tauks
#syn
    kout_one[4] = cellpar.alpha * FF * (1-ssc) - cellpar.beta * ssc

end #update_cell_I_br!

#This function calculates the steady-state variables for a specific V
#for FN neurons.
function steady_state_var_F(Varc_one, ion, cellpar, fl)
    Vc = Varc_one[1]

    hinf  = 1.0 / (1.0+exp((Vc-cellpar.thetah) / cellpar.sigmah)) 
    ninf  = 1.0 / (1.0+exp((Vc-cellpar.thetan) / cellpar.sigman))
    zinf  = 1.0 / (1.0+exp((Vc-cellpar.thetaz) / cellpar.sigmaz))
    uinf  = 1.0 / (1.0+exp((Vc-cellpar.thetau) / cellpar.sigmau))
    rinf  = 1.0 / (1.0+exp((Vc-cellpar.thetar) / cellpar.sigmar))

   
    FF = 1 / (1 + exp(-(Vc-cellpar.theta_syn)/cellpar.sigma_syn))

    Varc_one[2] = hinf
    Varc_one[3] = ninf
    Varc_one[4] = zinf
    Varc_one[5] = uinf
    Varc_one[6] = rinf
    Varc_one[7] = (cellpar.alpha * FF) / (cellpar.alpha * FF + cellpar.beta)
    
end #steady_state_var_F

#This function updates the variables of an IRT neuron model.
function update_cell_F!(Varc_one, kout_one, ion, cellpar, Isyn_cont, Iel_cont,
         s_ext, netpar, it, time, fl)
    Vc = Varc_one[1]
    hc = Varc_one[2]
    nc = Varc_one[3]
    zc = Varc_one[4]
    uc = Varc_one[5]
    rc = Varc_one[6]
    ssc = Varc_one[7]

    minf  = 1.0 / (1.0+exp((Vc-cellpar.thetam) / cellpar.sigmam))
    hinf  = 1.0 / (1.0+exp((Vc-cellpar.thetah) / cellpar.sigmah)) 
    tauh  = 30.0/(exp((Vc+50.0)/15.0)+exp(-(Vc+50.0)/16.0))
    INa   = cellpar.gNa  * minf*minf*minf * hc * (Vc-cellpar.VNa)

#    println(Vc, " ", hc, " ", nc, " ", zc, " ", uc, " ", rc , " ", ssc)
#    println(ion, " ", time, " ", INa, " ", minf, " ", nc, " ", Vc, " ", cellpar.VNa, " ",  cellpar.gNa)

    
    pinf  = 1.0 / (1.0+exp((Vc-cellpar.thetap) / cellpar.sigmap))
    INap  = cellpar.gNap * pinf * (Vc-cellpar.VNa)

    ninf  = 1.0 / (1.0+exp((Vc-cellpar.thetan) / cellpar.sigman))
    taun  = 7.0/(exp((Vc+40.0)/40.0)+exp(-(Vc+40.0)/50.0))
    IKdr  = cellpar.gKdr * nc*nc*nc*nc * (Vc-cellpar.VK)
    
    zinf  = 1.0 / (1.0+exp((Vc-cellpar.thetaz) / cellpar.sigmaz))
    IKm   = netpar.gKmar[ion] * zc * (Vc-cellpar.VK)

    uinf  = 1.0 / (1.0+exp((Vc-cellpar.thetau) / cellpar.sigmau))
    IAHP  = cellpar.gAHP * uc * (Vc-cellpar.VK)

    rinf  = 1.0 / (1.0+exp((Vc-cellpar.thetar) / cellpar.sigmar))
    taur  = 6000.0/(exp((Vc+140.0)/21.6)+exp(-(Vc+40.0)/22.7))
    Ih    = cellpar.gh * rc * (Vc-cellpar.Vh)

    IL   = netpar.gLar[ion] * (Vc-cellpar.VL)
  
    FF = 1 / (1 + exp(-(Vc-cellpar.theta_syn)/cellpar.sigma_syn))

#   Isyn_ext = netpar.gsyn_ext_ar[ion] * s_ext * (Vc - cellpar.Vsyn_ext)
    Isyn_ext = (netpar.gsyn_ext_ar[ion] * s_ext + cellpar.ginh_const) *
               (Vc - cellpar.Vsyn_ext) +
               cellpar.gexc_const * (Vc - cellpar.Vexc)

#V
    kout_one[1] = (-INa - INap - IKdr - IKm - IAHP - Ih -IL - Isyn_cont -
      Iel_cont - Isyn_ext + cellpar.Iapp) / cellpar.Cms

#h
    kout_one[2] = cellpar.phi*(hinf-hc)/tauh
#n
    kout_one[3] = cellpar.phi*(ninf-nc)/taun
#z
    kout_one[4] = (zinf-zc)/cellpar.tauz
#u
    kout_one[5] = (uinf-uc)/cellpar.tauu
#r
    kout_one[6] = (rinf-rc)/taur
#syn
    kout_one[7] = cellpar.alpha * FF * (1-ssc) - cellpar.beta * ssc

#   if ion==1 && ieq == 1
#       println("ion=", ion, " ieq=", ieq, " kout=", kout_one[1])
#   end
end #update_cell_F!

#This function calculates the steady-state variables for a specific V
#for I neurons based on the motorneuron (mn) model.
function steady_state_var_I_mn(Varc_one, ion, cellpar, fl)
    Vc = Varc_one[1]

    hinf  = 1.0 / (1.0+exp((Vc-cellpar.thetah) / cellpar.sigmah)) 
    ninf  = 1.0 / (1.0+exp((Vc-cellpar.thetan) / cellpar.sigman))
    zinf  = 1.0 / (1.0+exp((Vc-cellpar.thetaz) / cellpar.sigmaz))
    sinf  = 1.0 / (1.0+exp((Vc-cellpar.thetas) / cellpar.sigmas))
   
    FF = 1 / (1 + exp(-(Vc-cellpar.theta_syn)/cellpar.sigma_syn))

    Varc_one[2] = hinf
    Varc_one[3] = ninf
    Varc_one[4] = zinf
    Varc_one[5] = sinf
    Varc_one[6] = (cellpar.alpha * FF) / (cellpar.alpha * FF + cellpar.beta)
    
end #steady_state_var_I_mn

#This function updates the variables of an IRT neuron model.
function update_cell_I_mn!(Varc_one, kout_one, ion, cellpar, Isyn_cont,
         Iel_cont, s_ext, netpar, it, time, fl)
    Vc = Varc_one[1]
    hc = Varc_one[2]
    nc = Varc_one[3]
    zc = Varc_one[4]
    sc = Varc_one[5]
    ssc = Varc_one[6]

    minf  = Vc > -190.0 ? 1.0 / (1.0+exp((Vc-cellpar.thetam) / cellpar.sigmam)) : 0.0
    hinf  = Vc > -190.0 ? 1.0 / (1.0+exp((Vc-cellpar.thetah) / cellpar.sigmah)) : 1.0
    tauh  = 30.0/(exp((Vc+50.0)/15.0)+exp(-(Vc+50.0)/16.0))
    INa   = cellpar.gNa  * minf*minf*minf * hc * (Vc-cellpar.VNa)

#    println(Vc, " ", hc, " ", nc, " ", zc, " ", uc, " ", rc , " ", ssc)
#    println(ion, " ", time, " ", INa, " ", minf, " ", nc, " ", Vc, " ", cellpar.VNa, " ",  cellpar.gNa)
    
    pinf  = 1.0 / (1.0+exp((Vc-cellpar.thetap) / cellpar.sigmap))
    INap  = cellpar.gNap * pinf * (Vc-cellpar.VNa)

    ninf  = 1.0 / (1.0+exp((Vc-cellpar.thetan) / cellpar.sigman))
    taun  = 7.0/(exp((Vc+40.0)/40.0)+exp(-(Vc+40.0)/50.0))
    IKdr  = cellpar.gKdr * nc*nc*nc*nc * (Vc-cellpar.VK)
    
    zinf  = 1.0 / (1.0+exp((Vc-cellpar.thetaz) / cellpar.sigmaz))
    IKm   = netpar.gKmar[ion] * zc * (Vc-cellpar.VK)

    sinf  = 1.0 / (1.0+exp((Vc-cellpar.thetas) / cellpar.sigmas))
    ICas = cellpar.gCas * sc * (Vc - cellpar.VCa)

    IL   = netpar.gLar[ion] * (Vc-cellpar.VL)
  
    FF = 1 / (1 + exp(-(Vc-cellpar.theta_syn)/cellpar.sigma_syn))

#   Isyn_ext = netpar.gsyn_ext_ar[ion] * s_ext * (Vc - cellpar.Vsyn_ext)
    Isyn_ext = (netpar.gsyn_ext_ar[ion] * s_ext + cellpar.ginh_const) *
               (Vc - cellpar.Vsyn_ext) +
               cellpar.gexc_const * (Vc - cellpar.Vexc)
#   ion == 1 ? printfmtln("ion={1:d} gsyn_ext_ar={2:f} s_ext={3:f} Vc={4:f} Vsyn_ext={5:f} Isyn_ext={6:f}", ion, netpar.gsyn_ext_ar[ion], s_ext, Vc, cellpar.Vsyn_ext, Isyn_ext) : nothing

#V
    kout_one[1] = (-INa - INap - IKdr - IKm - ICas - IL - Isyn_cont -
      Iel_cont - Isyn_ext + cellpar.Iapp) / cellpar.Cms

#h
    kout_one[2] = cellpar.phi*(hinf-hc)/tauh
#n
    kout_one[3] = cellpar.phi*(ninf-nc)/taun
#z
    kout_one[4] = (zinf-zc)/cellpar.tauz
#s
    kout_one[5] = (sinf-sc)/cellpar.taus
#syn
    kout_one[6] = cellpar.alpha * FF * (1-ssc) - cellpar.beta * ssc

#   if ion==47
#       println(fl.trj, "ion=", ion, "V=", Vc, " time=", time, " zinf=", zinf,
#	" zc=", zc, " tauz=", cellpar.tauz)
#	println(fl.trj, "INa=", INa, " INap=", INap, " IKdr=", IKdr, " IKm=",
#	IKm, " ICas=", ICas, " IL=", IL, " Isyn_cont=", Isyn_cont, " Isyn_ext=",
#	Isyn_ext, " Iapp=", cellpar.Iapp)
#   end

end #update_cell_I_mn!

# This function substitues the functions related to each cell type in
# the corresponding cell structure.
function substitute_functions_in_cellpar(netpar)
    
    for istr in 1:length(netpar.C)
        expr = Meta.parse("netpar.C[" * string(istr) *
	"].steady_state_var = steady_state_var_" * netpar.C[istr].func_name)
#	println("istr=", istr, " expr=", expr)
        eval(expr)
        expr = Meta.parse("netpar.C[" * string(istr) *
	"].update_cell! = update_cell_" * netpar.C[istr].func_name * "!")
#	println("istr=", istr, " expr=", expr)
        eval(expr)
    end

end #substitute_functions_in_cellpar

# This function simulates a network with one parameter set.
function one_par(suffix, netpar, sval, av, partohop, rstpbstr, fl)

    netpar.C = deepcopy(netpar.C_orig)
    netpar.S = deepcopy(netpar.S_orig)
    
    [println("ipop=", ipop, " Iapp=", netpar.C[ipop].Iapp) for ipop in
      1:netpar.Npop]
    for ipop in 1:netpar.Npop, jpop in 1:netpar.Npop
        if netpar.connection_exists[ipop, jpop]
            println("ipop=", ipop, " jpop=", jpop, " gsyn=",
            netpar.S[ipop, jpop].syn_receptor_par_ar[1].gsyn)
	end
    end

    synstr = syn_str()
    statstr = stat_str()
    spkstr = spk_str()
    burstr = bur_str()
    
    construct_rst_pb_str(rstpbstr, netpar.Npop)
    
    runpar.tstat > runpar.Tall ? runpar.tstat = runpar.Tall : nothing
 
    substitute_functions_in_cellpar(netpar)

    substitute_values_in_structures(netpar, synstr, statstr, spkstr, runpar, av,
      fl)

    Varbar = Array{Array{Float64,1},1}(undef, netpar.Ntot)
    for ion in 1:netpar.Ntot
        ipop = netpar.whichpop[ion]
        Varbar[ion] = zeros(netpar.C[ipop].Neq)
    end

    #Substituting synaptic parameters of the connections to the post-synaptic
    #side
    initiate_synaptic_strengths_and_parameters(synstr, statstr,
      rstpbstr.prebot_burst_start, netpar, runpar, fl)

    #Substituting synaptic connectivity emerging from the pre-synaptic neuron.
    substitute_connectivity(netpar, runpar, fl)
    
    #Initial conditions
    in_con(Varbar, netpar, runpar, fl)
    
    #Integrating the differential euqations over the time interval
    n_run(Varbar, netpar, synstr, spkstr, statstr, runpar, rstpbstr, fl)   

    #Computing voltage statistics: chi
    compute_voltage_statistics(spkstr.Vav, statstr, netpar, runpar, av, fl)

    #Computing spike statistics: CV, CVt
    compute_spike_statistics(spkstr.ISIav, spkstr, netpar, runpar, av, fl)

    #Writing parameters on parthohop
    runpar.whisk_cal ? fill_in_par_to_hop(sval, netpar, runpar, partohop, fl) :
      nothing

    #Analazing bursting
    netpar.Npop >= 2 ? burst_analyze(burstr, rstpbstr.rst_ar, netpar.Ncells,
      runpar, av, fl) : nothing

end #one_par()

# This function constructs the arrays rst_ar and prebot_burst_start within
# the structure rstpbstr.
function construct_rst_pb_str(rstpbstr, Npop)

    rstpbstr.rst_ar = Array{Array{NamedTuple{(:tspk, :ion_pop),
      Tuple{Float64, Int64}},1}}(undef, Npop)
    [rstpbstr.rst_ar[ipop] = Array{NamedTuple{(:tspk, :ion_pop),
      Tuple{Float64, Int64}},1}(undef, 0) for ipop in 1:Npop]

    rstpbstr.prebot_burst_start = Array{Float64,1}(undef,0)

end #construct_rst_ar

# This function substitues values into the structures netpar, synstr, statstr
# and spkstr.
function substitute_values_in_structures(netpar, synstr, statstr, spkstr,
         runpar, av, fl)

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
    println("sm=", runpar.sm)
    
    multiply_K_N_by_Kfactor(netpar, fl)
    
    netpar.Ntot = sum(netpar.Ncells)
    netpar.pinds = zeros(Int64, netpar.Npop+1)
    netpar.whichpop = zeros(Int64, netpar.Ntot)

#   start index of each population (ends with netpar.Ntot+1)
    netpar.pinds = [1; ones(Int64, netpar.Npop) + cumsum(netpar.Ncells)]

    #simulation parameters

    netpar.whichpop = zeros(Int, netpar.Ntot)

    multiply_ginh_by_inhibition_factor(netpar, fl)
    
    netpar.gKmar = zeros(netpar.Ntot)
    netpar.gLar = zeros(netpar.Ntot)
    netpar.gsyn_ext_ar = zeros(netpar.Ntot)

#   vT = zeros(netpar.Ntot)
#   vTpop = [-40., -40., -40.]
#   sigvT = [3., 3., 3.]

    for ipop = 1:netpar.Npop
	non_min = netpar.pinds[ipop]
	non_max = netpar.pinds[ipop+1]-1

        onni = ones(Int64, netpar.Ncells[ipop])
        netpar.whichpop[non_min: non_max] =  ipop * onni

        onnr = ones(netpar.Ncells[ipop])

        netpar.gKmar[non_min: non_max] = netpar.C[ipop].gKm *
	onnr[1:netpar.Ncells[ipop]] + ((2.0 * rand(netpar.Ncells[ipop])) -
	onnr[1:netpar.Ncells[ipop]]) * netpar.C[ipop].del_gKm
		
        netpar.gLar[non_min: non_max] = netpar.C[ipop].gL *
	onnr[1:netpar.Ncells[ipop]] + ((2.0 * rand(netpar.Ncells[ipop])) -
	onnr[1:netpar.Ncells[ipop]]) * netpar.C[ipop].del_gL
		
	netpar.gsyn_ext_ar[non_min: non_max] = netpar.C[ipop].gsyn_ext *
	onnr[1:netpar.Ncells[ipop]] + ((2.0 * rand(netpar.Ncells[ipop])) -
	onnr[1:netpar.Ncells[ipop]]) * netpar.C[ipop].del_gsyn_ext

#       vT[non_min: non_max] = vTpop[ipop] * onnr[1:netpar.Ncells[ipop]] +
#       sigvT[ipop] * randn(netpar.Ncells[ipop])
    end

    if runpar.sm
        println(fl.out, "pinds=", netpar.pinds, " whichpop=", netpar.whichpop)
        for ion = 1:netpar.Ntot
            printfmtln(fl.out, "ion={1:d} gKmar={2:f} gLar={3:f} " *
	    "gsyn_ext_ar={4:f}", ion, netpar.gKmar[ion], netpar.gLar[ion],
	    netpar.gsyn_ext_ar[ion])
        end
    end	

    #Structure of synaptic variables
    synstr.savr_pop = zeros(netpar.Npop)
    synstr.Isyn_cont = zeros(netpar.Ntot)
    synstr.Iel_cont  = zeros(netpar.Ntot)
 
    #Structure of statistics of spikes, voltages and synaptic variables
    statstr.Vav_pop = zeros(netpar.Npop)
   
    #Structure of spikes.
#   spkstr.sp_in_dt = Array{Tuple{Float64, Int64, Int64}}(undef, netpar.Ntot)
#   The total number of spikes a neuron fired in a time step <= Ntot .
    spkstr.sp_in_dt = Array{NamedTuple{(:tspike, :jpop, :jon),
    Tuple{Float64, Int64, Int64}}}(undef, netpar.Ntot)

    spkstr.ndelay = zeros(netpar.Npop, netpar.Npop)
    spkstr.iptr_delay = zeros(netpar.Npop, netpar.Npop)
    spkstr.time_spike_delay = Array{Array{Array{Float64,1},1},1}(undef,
    netpar.Ntot)
    [spkstr.time_spike_delay[ion] = Array{Array{Float64,1},1}(undef,
    netpar.Npop) for ion=1:netpar.Ntot]

    spkstr.spk_during_prebot = zeros(Int64, netpar.Ntot)
    spkstr.spk_during_prebot_now = zeros(Int64, netpar.Ntot)
    spkstr.prebot_on_for_counting_spikes = false
    spkstr.n_prebot_on_for_counting_spikes = 0

    
    [netpar.connection_exists[ipop, jpop] ? compute_delay_values(
    netpar.S[ipop, jpop], spkstr, ipop, jpop, runpar, fl) : nothing for
    ipop=1:netpar.Npop for jpop=1:netpar.Npop]

    spkstr.Vav = Array{V_aver,1}(undef, netpar.Npop)
    [spkstr.Vav[ipop] = V_aver() for ipop=1:netpar.Npop]
    [spkstr.Vav[ipop].Vpop_avt = 0.0 for ipop=1:netpar.Npop]
    [spkstr.Vav[ipop].Vpop_sq_avt = 0.0 for ipop=1:netpar.Npop]
    [spkstr.Vav[ipop].V_avt = zeros(netpar.Ncells[ipop]) for ipop=1:netpar.Npop]
    [spkstr.Vav[ipop].V_sq_avt = zeros(netpar.Ncells[ipop])
      for ipop=1:netpar.Npop]
    [spkstr.Vav[ipop].t_spk_clamp = zeros(netpar.Ncells[ipop])
      for ipop=1:netpar.Npop]
    [spkstr.Vav[ipop].V_avt_spk_clamp = zeros(netpar.Ncells[ipop])
      for ipop=1:netpar.Npop]

    spkstr.ISIav = ISI_aver()
    for ion in 1:netpar.Ntot
        spkstr.ISIav.nspk = zeros(Int64, netpar.Ntot)
	spkstr.ISIav.tspk_old = zeros(netpar.Ntot)
	spkstr.ISIav.delta_tspk_old = zeros(netpar.Ntot)
        spkstr.ISIav.av  = zeros(netpar.Ntot)
        spkstr.ISIav.avt = zeros(netpar.Ntot)
        spkstr.ISIav.CV  = zeros(netpar.Ntot)
        spkstr.ISIav.CVt = zeros(netpar.Ntot)
    end
    for ipop in 1:netpar.Npop
        spkstr.ISIav.n_silent = zeros(Int64, netpar.Npop)
        spkstr.ISIav.CV_pop_av  = zeros(netpar.Npop)
        spkstr.ISIav.CVt_pop_av = zeros(netpar.Npop)
	spkstr.ISIav.T_pop_av = zeros(netpar.Npop)
    end
    
    for ipop = 1:netpar.Npop
        push!(av, avr_val())
    end
    av[1].Npop = netpar.Npop
    [av[ipop].Mfr = 0.0 for ipop in 1:netpar.Npop]
        
end #substitute_values_in_structure

# This function multiplies the inhibitory conductance strengths by
# inhibition_factor.
function multiply_ginh_by_inhibition_factor(netpar, fl)

    inhibition_threshold = -41.0

    println(fl.out, "inhibition_factor=", netpar.inhibition_factor)
    for ipop in 1:netpar.Npop
        if netpar.ar_cell_str[ipop] in netpar.inhibition_factor_consider &&
          netpar.C[ipop].Vsyn_ext < inhibition_threshold
            netpar.C[ipop].gsyn_ext = netpar.C_orig[ipop].gsyn_ext *
	      netpar.inhibition_factor
	    netpar.C[ipop].del_gsyn_ext = netpar.C_orig[ipop].del_gsyn_ext *
	      netpar.inhibition_factor
	end
    end

    println(fl.out, "\nmultiply by inhibition_factor")
    for ipop in 1:netpar.Npop
        printfmtln(fl.out, "ipop={1:d} Vsyn_ext={2:f} gsyn_ext={3:f} " *
	  "del_gsyn_ext={4:f}", ipop, netpar.C[ipop].Vsyn_ext,
	  netpar.C[ipop].gsyn_ext, netpar.C[ipop].del_gsyn_ext)
    end

    for ipop in 1:netpar.Npop
        for jpop in 1:netpar.Npop
	    if netpar.connection_exists[ipop, jpop]
	        if netpar.ar_cell_str[ipop] in 
                  netpar.inhibition_factor_consider &&
		  netpar.S[ipop, jpop].syn_receptor_par_ar[1].Vsyn <
		  inhibition_threshold
		    netpar.S[ipop, jpop].syn_receptor_par_ar[1].gsyn =
		      netpar.S_orig[ipop, jpop].syn_receptor_par_ar[1].gsyn *
		      netpar.inhibition_factor
		else
		    netpar.S[ipop, jpop].syn_receptor_par_ar[1].gsyn =
		      netpar.S_orig[ipop, jpop].syn_receptor_par_ar[1].gsyn
		end
            end
	end
    end

    for ipop in 1:netpar.Npop
        for jpop in 1:netpar.Npop
	    if netpar.connection_exists[ipop, jpop]
	        printfmtln(fl.out, "ipop={1:d} jpop={2:f} Vsyn={3:f} gsyn={4:f}",
	          ipop, jpop, netpar.S[ipop, jpop].syn_receptor_par_ar[1].Vsyn,
	          netpar.S[ipop, jpop].syn_receptor_par_ar[1].gsyn)
	    end
	end
    end

end #multiply_ginh_by_inhibition_factor

# This function multiplies the K's and the N's of the vIRT by Kfactor.
function multiply_K_N_by_Kfactor(netpar, fl)

    Kfactor_multiply = 2.0^netpar.Kfactor
    println(fl.out, "Kfactor=", netpar.Kfactor, " Kfactor_multiply=",
      Kfactor_multiply)
    
#   N
    netpar.Ncells = zeros(length(netpar.Ncells_orig))
    for ipop in 1:netpar.Npop
        if netpar.ar_cell_str[ipop][1:1] == "I"
	    netpar.Ncells[ipop] = trunc(Int64, netpar.Ncells_orig[ipop] *
	       Kfactor_multiply + runpar.epsilon)
	 else
	     netpar.Ncells[ipop] = netpar.Ncells_orig[ipop]
	end
    end

    println(fl.out, "\nmultiply by inhibition_factor")
    for ipop in 1:netpar.Npop
        printfmtln(fl.out, "ipop={1:d} Ncells={2:d}",
	  ipop, netpar.Ncells[ipop])
    end

#   K
    for ipop in 1:netpar.Npop
        for jpop in 1:netpar.Npop
	    if netpar.ar_cell_str[ipop][1:1] == "I" &&
	      netpar.ar_cell_str[jpop][1:1] == "I"
	        if netpar.connection_exists[ipop, jpop]
	            netpar.S[ipop, jpop].Kin = trunc(Int64,
		      netpar.S_orig[ipop, jpop].Kin *
		       Kfactor_multiply + runpar.epsilon)
	        end	  
            end
	end
    end
    
    for ipop in 1:netpar.Npop
        for jpop in 1:netpar.Npop
	    if netpar.connection_exists[ipop, jpop]
	        printfmtln(fl.out, "ipop={1:d} jpop={2:f} Kin={3:f}",
	          ipop, jpop, netpar.S[ipop, jpop].Kin)
	    end
	end
    end

end #multiply_K_N_by_Kfactor
    
#This function substitues delay values
function compute_delay_values(Sij, spkstr, ipop, jpop, runpar, fl)

#   Defining ndelay, the number of time steps in the delay stack.
    max_delay_time = Sij.tau_delay + Sij.Del_tau_delay
    ndelta_in_delay = trunc(Int64, (max_delay_time / runpar.deltat) +
    runpar.epsilon) + 2;
    spkstr.ndelay[ipop, jpop] = ndelta_in_delay;

#   Defining iptr_delay, the pointer indicating the present time step
#   in the stack.
    spkstr.iptr_delay[ipop, jpop] = spkstr.ndelay[ipop, jpop]

    [spkstr.time_spike_delay[ion][jpop] = zeros(spkstr.ndelay[ipop, jpop]) for
    ion in netpar.pinds[ipop]:(netpar.pinds[ipop+1]-1)]

    println(fl.out, "ipop=", ipop, " jpop=", jpop, " ndelay=",
    spkstr.ndelay[ipop, jpop], " iptr_delay=", spkstr.iptr_delay[ipop, jpop])

end #function compute_delay_value

#This function substitutes the connectivity matrix.
function substitute_connectivity(netpar, runpar, fl)

#   Chemical coupling
    netpar.nwcoup = Array{Array{NamedTuple{(:ipop, :ion_min, :ion_max),
            Tuple{Int64, Int64, Int64}}}}(undef, netpar.Ntot)
    netpar.wcoup = Array{Array{Int64,1}}(undef, netpar.Ntot)
    netpar.tdcoup = Array{Array{Float64,1}}(undef, netpar.Ntot)
    for jon in 1:netpar.Ntot
	netpar.nwcoup[jon] = []
        netpar.wcoup[jon] = []
        netpar.tdcoup[jon] = []
    end

    for ipop in 1:netpar.Npop
        for jpop in 1:netpar.Npop
            if netpar.connection_exists[ipop, jpop]
                find_coupling_matrix_zero_d(ipop, jpop, netpar, runpar, fl)
            end
        end
    end

#   Electrical coupling
   find_electrical_matrix_zero_d(netpar, runpar, fl)
    
end #substitute_connectivity

# This function finds the connectivity matrix for chemical synapses, sparse
# coupling and 0-d geometry.
function find_coupling_matrix_zero_d(ipop, jpop, netpar, runpar, fl)

    println(fl.out, "find_coupling_matrix_zero_d")
    println(fl.out, "ipop=", ipop, " range:", netpar.pinds[ipop], " ",
    netpar.pinds[ipop+1]-1)
    println(fl.out, "jpop=", jpop, " range:", netpar.pinds[jpop], " ",
    netpar.pinds[jpop+1]-1)

    for jon = netpar.pinds[jpop]:netpar.pinds[jpop+1]-1
        len_min = length(netpar.wcoup[jon]) + 1
        for ion = netpar.pinds[ipop]:netpar.pinds[ipop+1]-1
	    probx = netpar.S[ipop, jpop].Kin / netpar.Ncells[jpop]
	    if probx < 0.0 || probx > 1.0
	        println("probx=", probx, " is not between 0 and 1")
		exit(0)
	    end
	    rand_num = rand()
	    if rand_num < probx
	        push!(netpar.wcoup[jon], ion)
		
#               Finding tau_delay for this particular synaptic connection
		rand_num = rand()
		td_coup  = netpar.S[ipop, jpop].tau_delay +
	        netpar.S[ipop, jpop].Del_tau_delay * (2.0 * rand_num - 1.0)
	        push!(netpar.tdcoup[jon], td_coup)
	    end
	end
	len_max = length(netpar.wcoup[jon])

        if len_max - len_min >= 0
	    push!(netpar.nwcoup[jon], (ipop=ipop, ion_min=len_min,
	    ion_max=len_max))
	end
    end

    if runpar.sm
        println(fl.out, "ipop=", ipop, " jpop=", jpop)
        for jon = netpar.pinds[jpop]:netpar.pinds[jpop+1]-1
            println(fl.out, "jon=", jon, " len(wcoup)=", length(netpar.wcoup[jon]))
	    println(fl.out, "len(nwcoup)=", length(netpar.nwcoup[jon]))
	    println(fl.out, "nwcoup=", netpar.nwcoup[jon])
	    for icoup = 1:length(netpar.wcoup[jon])
	        printfmt(fl.out, "{1:d} {2:f} ", netpar.wcoup[jon][icoup],
	        netpar.tdcoup[jon][icoup])
	    end
	    println(fl.out, " ")
        end
    end
   
end #find_coupling_matrix_zero_d

# This function finds the connectivity matrix electrical coupling for
# 0-d geometry.
function find_electrical_matrix_zero_d(netpar, runpar, fl)
 
#   netpar.necoup = zeros(Int64, netpar.Ncells[ipop])
    netpar.ecoup = Array{Array{Tuple{Int64,Float64},1},1}(undef,netpar.Ntot)
    [netpar.ecoup[ion] = Array{Tuple{Int64,Float64}}(undef,0) for ion in 1:netpar.Ntot]
    
    for ion in 1:netpar.Ntot
        ipop = netpar.whichpop[ion]
	for jpop in ipop:netpar.Npop
	    if netpar.connection_exists[ipop, jpop]
                for jpre in 1:netpar.Ncells[jpop]
	            jon = netpar.pinds[jpop] -1 + jpre
		    if ion < jon &&
		          netpar.S[ipop,jpop].gel_one_psp > runpar.epsilon
	                probx = netpar.S[ipop, jpop].Kel / netpar.Ncells[jpop]
	                if probx < 0.0 || probx > 1.0
	                    println("probx=", probx, " is not between 0 and 1")
		        exit(0)
	                end
	                rand_num = rand()
	                if rand_num < probx
		            jon_gel_tuple = (jon,
			      netpar.S[ipop,jpop].gel_one_psp)
	      	            push!(netpar.ecoup[ion], jon_gel_tuple)
			    ion_gel_tuple = (ion,
			      netpar.S[ipop,jpop].gel_one_psp)
		            push!(netpar.ecoup[jon], ion_gel_tuple)
			end
		    end
                end
	    end
        end
    end

    if runpar.sm
        println(fl.out, "el:")
        for ion in 1:netpar.Ntot
	    print(fl.out, "ion=", ion, " necoup=", length(netpar.ecoup[ion]))
	    if length(netpar.ecoup[ion]) >= 1
	        print(fl.out, " ecoup=")
	        for jon_gel_tuple in netpar.ecoup[ion]
	            print(fl.out, jon_gel_tuple)
		    ch_pr = jon_gel_tuple != netpar.ecoup[ion][end] ?
		      " " : "\n"
		    print(fl.out, ch_pr)
	        end
	    else
	        print(fl.out, "\n")
	    end
	end
	println(fl.out, "  ")
	
	for ipop in 1:netpar.Npop
	    ne_av = 0.0
	    for ion in netpar.pinds[ipop]:netpar.pinds[ipop+1]-1
	        ne_av += length(netpar.ecoup[ion])
	    end
	    ne_av /= netpar.Ncells[ipop]
	    println(fl.out, "ipop=", ipop, " ne_av=", ne_av)
	end
    end

end #find_electrical_matrix_zero_d


# This function substitues the values of the synaptic synaptic parameters
# of the synapses to the post-synaptic neurons into the structure synstr.
function initiate_synaptic_strengths_and_parameters(synstr, statstr,
  prebot_burst_start, netpar, runpar, fl)

    synstr.nsynvar = zeros(Int64, netpar.Npop)
    synstr.synpar_during_run = Array{Array{syn_par_during_run,1},1}(
    undef,netpar.Npop)
    synstr.isyn_to_send = Array{Array{Array{Int64,1},1},1}(undef,netpar.Npop)

#   Chemical synapses
    for ipop=1:netpar.Npop
	synstr.synpar_during_run[ipop] = Array{syn_par_during_run,1}(undef, 0)
	synstr.isyn_to_send[ipop] = Array{Array{Int64,1},1}(undef,netpar.Npop)

        for jpop=1:netpar.Npop
	    if netpar.connection_exists[ipop, jpop]
                synstr.isyn_to_send[ipop][jpop] = zeros(Int64, 2)

	        isyn_to_send_1 = synstr.nsynvar[ipop] + 1
	        len_ar = length(netpar.S[ipop, jpop].syn_receptor_par_ar)
	        for ireceptor in 1:len_ar
#                   This function also defines nsynvar. 		
	            one_synaptic_type_strengths_parameters(netpar.S[ipop,
		    jpop], ipop, jpop, synstr, runpar, fl)
		end
		isyn_to_send_2 = synstr.nsynvar[ipop]

                if isyn_to_send_2 >= isyn_to_send_1
		    synstr.isyn_to_send[ipop][jpop][1] = isyn_to_send_1
		    synstr.isyn_to_send[ipop][jpop][2] = isyn_to_send_2
		else
		    println("ion=", ion, " jon=", jon,
		    " isyn_to_send_2 < isyn_to_send_1")
		    exit(0)
		end
	    end
	end
    end

#   Electrical coupling
    for ipop=1:netpar.Npop
        for jpop=1:netpar.Npop
	    println(fl.out, "ipop=", ipop, " jpop=", jpop)
	    if netpar.connection_exists[ipop, jpop]
                netpar.S[ipop,jpop].gel_one_psp = netpar.S[ipop,jpop].gel /
                  netpar.Kel_scale
	    end
	end
    end


    if runpar.sm
        println(fl.out, "\nsynpar_during_run")
        for ipop=1:netpar.Npop
            println(fl.out, "ipop=", ipop, " nsynvar=", synstr.nsynvar[ipop])
            for isynvar=1:synstr.nsynvar[ipop]
	        println(fl.out, ipop, " ", isynvar, " ",
	        synstr.synpar_during_run[ipop][isynvar].gsyn_one_psp, " ",
	        synstr.synpar_during_run[ipop][isynvar].Vsyn, " ",
	        synstr.synpar_during_run[ipop][isynvar].tsyn, " ",
	        synstr.synpar_during_run[ipop][isynvar].is_it_NMDA)
	    end
        end

        println(fl.out, "isyn_to_send")
        for ipop=1:netpar.Npop, jpop=1:netpar.Npop
	    netpar.connection_exists[ipop, jpop] ? println(fl.out, ipop, " ",
	    jpop, " ", synstr.isyn_to_send[ipop][jpop][1]) : nothing
        end

        println(fl.out, "gel_one_psp_during_run")
        for ipop=1:netpar.Npop
            for jpop=1:netpar.Npop
 	        if netpar.connection_exists[ipop, jpop]
                    println(fl.out, ipop, " ", jpop, " ",
		      netpar.S[ipop,jpop].gel_one_psp)
		end
	    end
	end
        println(fl.out, " ")
    end
    
#   Defining synvar
    synstr.synvar = Array{Array{Float64,1},1}(undef,netpar.Ntot)
    for ion in 1:netpar.Ntot
        ipop = netpar.whichpop[ion]
        synstr.synvar[ion] = zeros(synstr.nsynvar[ipop])
    end

    #Structure of statistics of spikes, voltages and synaptic variables
    statstr.sav_pop = Array{Array{Float64,1},1}(undef,netpar.Npop)
    if runpar.consider_s == 'r'
        [statstr.sav_pop[ipop] = zeros(1) for ipop=1:netpar.Npop]
    elseif runpar.consider_s == 's'
        [statstr.sav_pop[ipop] = zeros(synstr.nsynvar[ipop]) for
	  ipop=1:netpar.Npop]
    end

    statstr.Sp_av  = 0.0
    statstr.Sp_avt = 0.0
    statstr.Sm_av  = 0.0
    statstr.Sm_avt = 0.0
    
#   Initiates variables for external synaptic strengths.
    synstr.s_ext = zeros(netpar.Npop)
    synstr.s_ext_previous = zeros(netpar.Npop)
    compute_random_Tper(synstr, 0.0, 0, prebot_burst_start, netpar, runpar, fl)

end #initiate_synaptic_strengths_and_parameters

#This function computes the next value of Tper_syn_ext_rand for
#random type_syn_ext.
function compute_random_Tper(synstr, time, it, prebot_burst_start, netpar,
  runpar, fl)

    if netpar.type_syn_ext == 'p'
        synstr.Tper_syn_ext_rand = netpar.Tper_syn_ext
    elseif netpar.type_syn_ext == 'r'
        rand_num = rand()
        Tper_rand = netpar.Tper_syn_ext +
	  0.5 * netpar.Trand_syn_ext * (2.0 * rand_num - 1.0)
	synstr.Tper_syn_ext_rand = floor(Tper_rand / runpar.deltat) *
	  runpar.deltat
    else
        println("type_syn_ext=", netpar.type_syn_ext, " should be p or r")
	exit(0)
    end

    synstr.T_start_cycle = time

    push!(prebot_burst_start, time)
    if runpar.sm
        printfmtln(fl.out, "it={1:d} time={2:f} Tper_syn_ext_rand={3:f}" *
        " T_start_cycle={4:f}", it, time, synstr.Tper_syn_ext_rand, 
        synstr.T_start_cycle)
    end
      
end #compute_random_Tper

# This function initiates the structure of synaptic variables for one
# type of synapses (AMPA, NMDA, GABAA).
function one_synaptic_type_strengths_parameters(Sij, ipop, jpop, synstr, runpar,
         fl)

    synstr.nsynvar[ipop] += 1
#   synpar_during_run = syn_par_during_run()
    push!(synstr.synpar_during_run[ipop], syn_par_during_run())
    
    for ion in netpar.pinds[ipop]:netpar.pinds[ipop+1]-1
        synstr.synpar_during_run[ipop][synstr.nsynvar[ipop]].Vsyn =
	  Sij.syn_receptor_par_ar[1].Vsyn
        synstr.synpar_during_run[ipop][synstr.nsynvar[ipop]].tsyn =
	  Sij.syn_receptor_par_ar[1].tsynd
        synstr.synpar_during_run[ipop][synstr.nsynvar[ipop]].is_it_NMDA =
	  Sij.syn_receptor_par_ar[1].is_it_NMDA
	if !synstr.synpar_during_run[ipop][synstr.nsynvar[ipop]].is_it_NMDA
	  synstr.synpar_during_run[ipop][synstr.nsynvar[ipop]].gsyn_one_psp =
	  Sij.syn_receptor_par_ar[1].gsyn /
	  (netpar.S[ipop, jpop].Kin * Sij.syn_receptor_par_ar[1].tsynd)
	else
	  # correct this!!!
	  synstr.synpar_during_run[ipop][synstr.nsynvar[ipop]].gsyn_one_psp =
	  Sij.syn_receptor_par_ar[1].gsyn /
	  (netpar.S[ipop, jpop].Kin * netpar.S[ipop, jpop].tsynd)
            synstr.synpar_during_run[ipop][synstr.nsynvar[ipop]].thetanp =
  	      Sij.syn_receptor_par_ar[1].thetanp
            synstr.synpar_during_run[ipop][synstr.nsynvar[ipop]].sigmanp =
	      Sij.syn_receptor_par_ar[1].sigmanp
        end
    end

    if Sij.syn_receptor_par_ar[1].is_it_NMDA
        synstr.nsynvar[ipop] += 1
	push!(synstr.synpar_during_run[ipop],
	  synstr.synpar_during_run[ipop][synstr.nsynvar[ipop]-1])
	synstr.synpar_during_run[ipop][synstr.nsynvar[ipop]].tsyn =
	  Sij.syn_receptor_par_ar[1].tsynr
    end

#   computing the extremum of single PSPs using a linear approximation.    
    NMDA_is = 0
    if Sij.syn_receptor_par_ar[1].is_it_NMDA
        NMDA_is = 1
        ftau = functau_NMDA(netpar.C[ipop].gL, netpar.C[ipop].Cms,
	  Sij.syn_receptor_par_ar[1].tsynr,
	  Sij.syn_receptor_par_ar[1].tsynd, runpar)
    else
        ftau = functau(netpar.C[ipop].gL, netpar.C[ipop].Cms,
	Sij.syn_receptor_par_ar[1].tsynd)
    end
    Vextr = -(
      synstr.synpar_during_run[ipop][synstr.nsynvar[ipop]].gsyn_one_psp *
      netpar.S[ipop, jpop].UU *
        (netpar.C[ipop].VL -
	synstr.synpar_during_run[ipop][synstr.nsynvar[ipop]].Vsyn) /
	netpar.C[ipop].Cms) * ftau * (Sij.syn_receptor_par_ar[1].tsynd -
	Sij.syn_receptor_par_ar[1].tsynr)

    if runpar.sm
        printfmtln(fl.out, "ipop={1:d} jpop={2:d} gsyn_one_psp={3:f} " *
	    "ftau={4:f} NMDA_is={5:d} Cm={6:f} Vextr={7:f}", ipop, jpop,
            synstr.synpar_during_run[ipop][synstr.nsynvar[ipop]].gsyn_one_psp,
            ftau, NMDA_is, netpar.C[ipop].Cms, Vextr)
    end
    
#   printfmtln("UU={1:f} VL={2:f} Vsyn={3:f} tsynd={4:f} tsynr={5:f}",
#       netpar.S[ipop, jpop].UU, netpar.C[ipop].VL,
#       synstr.synpar_during_run[ipop][synstr.nsynvar[ipop]].Vsyn,
#	Sij.syn_receptor_par_ar[1].tsynd, Sij.syn_receptor_par_ar[1].tsynr)

end #one_synaptic_type_strengths_parameters

#This function compute the function of tau0=Cm/gL and taus=tsyn that
#yields the exremum PSP for an exponential PSC.
function functau(gL, Cm, tsyn)
  tau0 = Cm / gL
  taus = tsyn
  dft = tau0 - taus
  rtt = taus / tau0
  tp = -tau0 * taus * log(rtt) / dft
  ftau = (tau0 / dft) * (rtt^(taus / dft) - rtt^(tau0 / dft))

  return ftau
end #functau

#This function compute the maximim of the function of tau0=Cm/gL tausa
#and tausb for extremem NMDA PSP.
function functau_NMDA(gL, Cm, tausa, tausb, runpar)
    numt = 1000
    tau0 = Cm / gL
    tmax = 10.0 * max(tau0, tausa, tausb)

    if abs(tau0 - tausa) < runpar.epsilon
        printfmtln("tau0={1:f} = tausa={2:f}", tau0, tausa)
        exit(0)
    elseif abs(tau0 - tausb) < runpar.epsilon
        printfmtln("tau0={1:f} = tausb={2:f}", tau0, tausb)
        exit(0)
    elseif abs(tausa - tausb) < runpar.epsilon
        printfmtln("tausa={1:f} = tausb={2:f}", tausa, tausb)
        exit(0);
    end
  
    coef0 = tau0 * tau0  / ((tausa - tau0) * (tausb - tau0))
    coef1 = tau0 * tausa / ((tausa - tau0) * (tausb - tausa))
    coef2 = tau0 * tausb / ((tausb - tau0) * (tausb - tausa))

    itmax = 0;
    Vfun_max = 0.0;

    for it in 1:numt
        tt = it * tmax / numt;
        Vfun = coef0 * exp(-tt / tau0) - coef1 * exp(-tt / tausa) + 
          coef2 * exp(-tt / tausb)
        if (Vfun > Vfun_max)
            itmax = it
            Vfun_max = Vfun
        end
    end

    if (itmax >= 2) && (itmax < numt-1)
        t0 = (itmax-1) * tmax / numt
        t1 = itmax * tmax / numt
        t2 = (itmax+1) * tmax / numt
        V0 = coef0 * exp(-t0 / tau0) - coef1 * exp(-t0 / tausa) + 
          coef2 * exp(-t0 / tausb)
        V1 = coef0 * exp(-t1 / tau0) - coef1 * exp(-t1 / tausa) + 
          coef2 * exp(-t1 / tausb)
        V2 = coef0 * exp(-t2 / tau0) - coef1 * exp(-t2 / tausa) + 
          coef2 * exp(-t2 / tausb)

        xb = V2 - V0
        xc = V0 - 2.0 * V1 + V2
        if abs(xc) < runpar.epsilon
            Vfun_max = V1
        else        
           Vfun_max = V1 - 0.125 * xb * xb / xc
        end
    end

    ftau  = Vfun_max
    return ftau
end #functau_NMDA

#This function substitutes the initial conditions.
function in_con(Varbar, netpar, runpar, fl)

    if (runpar.incond == 'a') # for checking
       for ion = 1:netpar.Ntot
            ipop = netpar.whichpop[ion]
            Varbar[ion][1] = netpar.C[ipop].Vinc1 + (netpar.C[ipop].Vinc2 -
	      netpar.C[ipop].Vinc1) * ion * 0.01
       end
    elseif (runpar.incond == 'r') # read initial conditions
        for ion = 1:netpar.Ntot
	    ipop = netpar.whichpop[ion]
	    [Varbar[ion][ieq] = netpar.C[ipop].initial_val[ieq] for
	      ieq=1:netpar.C[ipop].Neq]
        end    
    elseif (runpar.incond == 'b') # auxiliary variables calculated from V
        for ion = 1:netpar.Ntot
	    ipop = netpar.whichpop[ion]
            Varbar[ion][1] = netpar.C[ipop].Vinc1 +
	    (netpar.C[ipop].Vinc2 - netpar.C[ipop].Vinc1) * rand()
            netpar.C[ipop].steady_state_var(Varbar[ion], ion, netpar.C[ipop],
	      fl)
        end
    elseif (runpar.incond == 's')  # steady state
        for ion = 1:netpar.Ntot
	    ipop = netpar.whichpop[ion]
            Varbar[ion][1] = netpar.C[ipop].Vinc1
            netpar.C[ipop].steady_state_var(Varbar[ion], ion, netpar.C[ipop],
	      fl)
	    Varb = Varbar[ion]
	    update_cell_modified!(KK, VV) = netpar.C[ipop].update_cell!(VV,
	      KK, ion, netpar.C[ipop], 0.0, 0.0, netpar, 0, 0.0, fl)
	    sol_steady_state = nlsolve(update_cell_modified!, Varb)
#	    println(sol_steady_state.zero)
            Varbar[ion] = sol_steady_state.zero
	    # To escape from an unstable FP.
	    Varbar[ion][1] += runpar.Vincond_add + runpar.Vincond_rand * rand()
        end
    end
end #in_con

function n_run(Varbar, netpar, synstr, spkstr, statstr, runpar, rstpbstr, fl)

    k0 = Array{Array{Float64,1},1}(undef, netpar.Ntot)
    [k0[ion] = zeros(netpar.C[netpar.whichpop[ion]].Neq) for ion=1:netpar.Ntot]

    k1 = Array{Array{Float64,1},1}(undef, netpar.Ntot)
    [k1[ion] = zeros(netpar.C[netpar.whichpop[ion]].Neq) for ion=1:netpar.Ntot]
    
    k2 = Array{Array{Float64,1},1}(undef, netpar.Ntot)
    [k2[ion] = zeros(netpar.C[netpar.whichpop[ion]].Neq) for ion=1:netpar.Ntot]
    
    k3 = Array{Array{Float64,1},1}(undef, netpar.Ntot)
    [k3[ion] = zeros(netpar.C[netpar.whichpop[ion]].Neq) for ion=1:netpar.Ntot]
    
    k4 = Array{Array{Float64,1},1}(undef, netpar.Ntot)
    [k4[ion] = zeros(netpar.C[netpar.whichpop[ion]].Neq) for ion=1:netpar.Ntot]
    
    Varc = Array{Array{Float64,1},1}(undef, netpar.Ntot)
    [Varc[ion] = zeros(netpar.C[netpar.whichpop[ion]].Neq) for ion=1:netpar.Ntot]
    
    Varold = zeros(netpar.Ntot, 2)
    after_max_vol = zeros(Int64, netpar.Ntot)

    runpar.first_print = true
    
    #simulation-------------------------------------------------------------
    println("starting sim")
    println("NT=", runpar.NT)

    it = 0
    time = it * runpar.deltat

    compute_statistics(Varbar, synstr, statstr, netpar, runpar, fl)
    
    if runpar.sm && runpar.tmcol + runpar.epsilon >= runpar.Tall
        pr_fct(Varbar, statstr, synstr, netpar, runpar, time, it, fl)
    end

    for it = 1:runpar.NT
        time = it * runpar.deltat
	if Int(mod(it, runpar.NT/2)) == 1 #print percent complete
	    println(round(Int,100*it/ runpar.NT),"%")
	    #print("\r",round(Int,100*it/ runpar.NT),"%")
	end

	for ion_all = 1:netpar.Ntot #update synaptic/adaptation parameters
	    ion = netpar.whichpop[ion_all]
	end

	#Computing chemical synaptic currents
	if runpar.consider_s == 'r'
	    compute_total_synaptic_conductance_on_a_neuron_pre(Varbar, synstr,
	      time, it, netpar, runpar, fl)
	elseif runpar.consider_s == 's'
	    compute_total_synaptic_conductance_on_a_neuron_post(Varbar, synstr,
	      time, it, netpar, runpar, fl)
	else
	    println("consider_s should be r or s")
	    exit(0)
	end

        #Computing electrical synaptic currents
	compute_total_electric_conductance_on_a_neuron(Varbar, synstr, time, it,
	  netpar, runpar, fl)

        #Computing external synaptic variable
	external_gsyn_value(synstr, time, it, rstpbstr.prebot_burst_start,
	  netpar, runpar, fl)

        #Advancing the state vetor for one time step.
	if runpar.method  == 'r'                #Runge-Kutta-4 method
            one_integration_step(Varbar, k0, k1, Varc, synstr, 0.0            ,
	      time, it, netpar.Ntot, netpar, runpar, fl)
            one_integration_step(Varbar, k1, k2, Varc, synstr, runpar.deltat/2,
	      time, it, netpar.Ntot, netpar, runpar, fl)
            one_integration_step(Varbar, k2, k3, Varc, synstr, runpar.deltat/2,
	      time, it, netpar.Ntot, netpar, runpar, fl)
            one_integration_step(Varbar, k3, k4, Varc, synstr, runpar.deltat  ,
	      time, it, netpar.Ntot, netpar, runpar, fl)

            for ion = 1:netpar.Ntot
	        for ieq = 1:netpar.C[netpar.whichpop[ion]].Neq 
	            Varbar[ion][ieq] += (runpar.deltat/6.0) * (k1[ion][ieq] +
		      2 * k2[ion][ieq] + 2 * k3[ion][ieq] + k4[ion][ieq])
	        end
            end
	elseif runpar.method  == 't'            #Runge-Kutta-2 method
            one_integration_step(Varbar, k0, k1, Varc, synstr, 0.0            ,
	      time, it, netpar.Ntot, netpar, runpar, fl)
            one_integration_step(Varbar, k1, k2, Varc, synstr, runpar.deltat/2,
	      time, it, netpar.Ntot, netpar, runpar, fl)
	    
            for ion = 1:netpar.Ntot
	        for ieq = 1:netpar.C[netpar.whichpop[ion]].Neq
	            Varbar[ion][ieq] += runpar.deltat * k2[ion][ieq]
	        end
            end
	elseif runpar.method  == 'e'            #Euler method 
            one_integration_step(Varbar, k0, k1, Varc, synstr, 0.0            ,
	      time, it, netpar.Ntot, netpar, runpar, fl)
	    
	    for ion = 1:netpar.Ntot
	        for ieq = 1:netpar.C[netpar.whichpop[ion]].Neq
		    Varbar[ion][ieq] += runpar.deltat * k1[ion][ieq]
		end
	    end
	end

        #Checking for nan.
	for ion in 1:netpar.Ntot
	    if any(isnan, Varbar[ion])
	        println("NaN in Varbar[ion]!", ion)
	        exit(0)
	    end
	end

	#Noise
	if netpar.noise >= runpar.epsilon
	    Varbar[:][1] += sqrt(2.0 * netpar.noise * runpar.deltat) * randn(netpar.Ntot)
	end

        #Spike detection
	if runpar.generate_one_psp == 'n' #Network simulation
	    spike_detect(Varbar, Varold, after_max_vol, it, time,
	      rstpbstr.rst_ar, synstr, spkstr, netpar, runpar, fl)
	    if time >= runpar.Tall - runpar.tstat + 0.5 * runpar.deltat
	        spikes_during_prebot_activity(it, time, synstr, spkstr, netpar,
	          runpar, fl)
            end
	elseif runpar.generate_one_psp == 'y' #Generate only single PSPs
	    spike_generate(Varbar, Varold, after_max_vol, it, time,
	    rstpbstr.rst_ar, synstr, spkstr, netpar, runpar, fl)
        else
	    println("wrong generate_one_psp=", runpar.generate_one_psp)
	    exit(0)
	end

        #Spread of spike timing to post-synaptic neurons with time delay.
        multiple_store_spikes_plus_td(it, time, synstr, spkstr, netpar, runpar,
	  fl)

        #Post-synaptic variables s decay exponentially
        decay_post_synaptic_variables(synstr, time, it, netpar, runpar, fl)

        update_post_synaptic_variables_for_pre_synaptic_spikes(synstr, spkstr,
	  time, it, netpar, runpar, fl)
	
	for ion = 1:netpar.Ntot
            Varold[ion, 2] = Varold[ion, 1]
            Varold[ion, 1] = Varbar[ion][1]
	end
        [synstr.s_ext_previous[ipop] = synstr.s_ext[ipop] for ipop in
	  1:netpar.Npop]
	
	if time >= runpar.Tall - runpar.tstat + 0.5 * runpar.deltat
	    update_Vav_arrays(Varbar, spkstr.Vav, it, time, runpar.deltat,
	    netpar, runpar, fl)
	    update_ISIav_arrays(spkstr, it, time, netpar, runpar, fl)
	end
	
        compute_statistics(Varbar, synstr, statstr, netpar, runpar, fl)

        no_sav_pop = length(statstr.sav_pop) == 1 &&
	  length(statstr.sav_pop[1]) == 0
	if !no_sav_pop
	    if time >= runpar.Tall - runpar.tstat + 0.5 * runpar.deltat
	        compute_syn_statistics(statstr, netpar, runpar, fl)
	    end
        end
	
	if runpar.sm && time >= runpar.Tall -runpar.tmcol + runpar.epsilon &&
	  it%runpar.twrite == 0
            pr_fct(Varbar, statstr, synstr, netpar, runpar, time, it, fl)
	end

    end
end #n_run

#This function computes one integration step.
function one_integration_step(Varbar, kin, kout, Varc, synstr, delt, time, it,
         Ntot, netpar, runpar, fl)

#   Runge-Kutta input variables
    if (delt > runpar.epsilon)
        Varc = Varbar + delt * kin
    else
        Varc = Varbar
    end

    for ion = 1:netpar.Ntot
        ipop = netpar.whichpop[ion]
	netpar.C[ipop].update_cell!(Varc[ion], kout[ion], ion,
	  netpar.C[ipop], synstr.Isyn_cont[ion], synstr.Iel_cont[ion],
	  synstr.s_ext[ipop], netpar, it, time, fl)
    end
    
#   println(fl.col, time, " Isyn1=", synstr.Isyn_cont[1])
end #one_integration_step()

#This function computes the chemical synaptic fields on the neurons.
#The synaptic variables are pre-synaptic.
function compute_total_synaptic_conductance_on_a_neuron_pre(Varbar, synstr,
         time, it, netpar, runpar, fl)

    #All-to-all
    for ipop = 1:netpar.Npop
        isvar = netpar.C[ipop].Neq
	synstr.savr_pop[ipop] = 0.0
	for ion in netpar.pinds[ipop]:netpar.pinds[ipop+1]-1
	    synstr.savr_pop[ipop] += Varbar[ion][isvar]
	end
	synstr.savr_pop[ipop] /= netpar.pinds[ipop+1] - netpar.pinds[ipop]
##      synstr.savr_pop[ipop] =
##        mean(Varbar[netpar.pinds[ipop]:(netpar.pinds[ipop+1]-1)][isvar])
#	if it == 1000
#   	    println(fl.out, "it=", it, " ipop=", ipop, " savr_pop=",
#             synstr.savr_pop[ipop])
#	end
    end

    for ion = 1:netpar.Ntot
        ipop = netpar.whichpop[ion]
        synstr.Isyn_cont[ion] = 0.0
	for jpop = 1:netpar.Npop
            if netpar.connection_exists[ipop, jpop]
	        gsyn_one_psp =
		  netpar.S[ipop, jpop].syn_receptor_par_ar[1].gsyn /
		  netpar.S[ipop, jpop].Kin
		Vsyn = netpar.S[ipop, jpop].syn_receptor_par_ar[1].Vsyn
		DelV = netpar.rho_concur * (Varbar[ion][1] - Vsyn) +
		  (1.0 - netpar.rho_concur) *  (netpar.C[ipop].VL - Vsyn)
                Isyn_from_one_pop = gsyn_one_psp *  synstr.savr_pop[jpop] *
	          DelV
#  	        if it == 50100 && ion == 1
#	            println(fl.out, "ion=", ion, " ipop=", ipop, " jpop=", jpop,
#                   " gsyn=", gsyn, " savr=", synstr.savr_pop[jpop], " V=",
#	  	    Varbar[ion][1], " ", Vsyn, " Isyn=", Isyn_from_one_pop)
#               end
                synstr.Isyn_cont[ion] += Isyn_from_one_pop
            end
        end

        if it == 1000
            println(fl.out, "ion=", ion, " Isyn_cont=", synstr.Isyn_cont[ion])
        end
    end
end #compute_total_synaptic_conductance_on_a_neuron_pre

#This function computes the chemical synaptic fields on the neurons.
#The synaptic variables are post-synaptic.
function compute_total_synaptic_conductance_on_a_neuron_post(Varbar, synstr,
         time, it, netpar, runpar, fl)

#   println(time, " ", synstr.synvar[1][1])
#   Chemical synapses
    for ion in 1:netpar.Ntot
        ipop = netpar.whichpop[ion]
        Vpost = Varbar[ion][1]
	synstr.Isyn_cont[ion] = 0.0

        for isynvar in 1:synstr.nsynvar[ipop]
	   synstr_here = synstr.synpar_during_run[ipop][isynvar]
	    fNMDA = !synstr_here.is_it_NMDA ?
	      1.0 : Gammaf(Vpost, synstr_here.thetanp, synstr_here.sigmanp)
	    DelV = netpar.rho_concur * (Vpost - synstr_here.Vsyn) +
	      (1.0 - netpar.rho_concur) * (netpar.C[ipop].VL - synstr_here.Vsyn)
#	    if ion == 1
#	        println(it, " ", time, " ", Vpost, " ", " ion=", ion, " ipop=", ipop, " ns=", synstr.nsynvar[ipop], " gs=", synstr_here.gsyn_one_psp, " DelV=", DelV, " fNMDA=", fNMDA, " s=", synstr.synvar[ion][isynvar])
#	    end
	    synstr.Isyn_cont[ion] += synstr_here.gsyn_one_psp *
	      synstr.synvar[ion][isynvar] * DelV * fNMDA
	end

#  	if it == 50100 && ion == 1
#	    println(fl.out, "ion=", ion, " ipop=", ipop, 
#           " gsyn=", gsyn, " savr=", synstr.savr_pop[jpop], " V=",
#	    Varbar[ion][1], " ", Vsyn, " Isyn=", Isyn_from_one_pop)
#        end
    end

#    println(time, " ", synstr.Isyn_cont[1], " ", synstr.synvar[1][1], " ", synstr.synvar[1][2])
end #compute_total_synaptic_conductance_on_a_neuron_post

#This function computes the electrical synaptic fields on the neurons.
function compute_total_electric_conductance_on_a_neuron(Varbar, synstr,time, it,
         netpar, runpar, fl)

    for ion in 1:netpar.Ntot
        synstr.Iel_cont[ion] = 0.0
	for e_tuple in netpar.ecoup[ion]
	    jon = e_tuple[1]
	    gel_one_psp = e_tuple[2]
	    Vdiff = Varbar[ion][1] - Varbar[jon][1]
	    synstr.Iel_cont[ion] += gel_one_psp * Vdiff
	end
    end

end #compute_total_electric_conductance_on_a_neuron

#This function compute the values of the external synaptic variables.
function external_gsyn_value(synstr, time, it, prebot_burst_start, netpar,
  runpar, fl)
   
    t_in_ext_syn_cycle = time - runpar.deltat + runpar.epsilon -
      synstr.T_start_cycle
    s_ext_all = t_in_ext_syn_cycle < netpar.Tup_syn_ext ? 1.0 : 0.0
    
    for ipop in 1:netpar.Npop
        synstr.s_ext[ipop] = s_ext_all
    end

    if runpar.sm
        ipop = 1
	if it == 1
	    printfmtln(fl.sex, "{1:f} {2:f}", time - runpar.deltat,
              synstr.s_ext[ipop])
	elseif abs(synstr.s_ext[ipop] - synstr.s_ext_previous[ipop]) >
	  runpar.epsilon
	    printfmtln(fl.sex, "{1:f} {2:f}", time - runpar.deltat,
	      synstr.s_ext_previous[ipop])
            printfmtln(fl.sex, "{1:f} {2:f}", time - runpar.deltat,
	      synstr.s_ext[ipop])
	end
    end
    
    t_for_cal_new_T_start = time + runpar.epsilon - synstr.T_start_cycle
    if t_for_cal_new_T_start > synstr.Tper_syn_ext_rand
        compute_random_Tper(synstr, time, it, prebot_burst_start, netpar,
	  runpar, fl)
    end
    
end #external_gsyn_value

#This function finds the spikes fired by pre-synaptic neurons and updates
#the synaptic variables of the post-synaptic neurons.
function spike_detect(Varbar, Varold, after_max_vol, it, time, rst_ar, synstr,
         spkstr, netpar, runpar, fl)

    tpeak = -1.0
    Vpeak = -100.0
    spkstr.nsp_in_dt = 0
    
    for jon = 1:netpar.Ntot
        V0 = Varbar[jon][1]
	if after_max_vol[jon] == 0 && it >= 10 
	    V1 = Varold[jon][1]

            if runpar.spike_detect_criterion == 'p'
                V2 = Varold[jon][1]
                spike_detected, tpeak, Vpeak = spike_detect_peak(V0, V1,
		V2, time, tpeak, Vpeak, netpar, runpar, fl)
            elseif runpar.spike_detect_criterion == 't'
                spike_detected, tpeak, Vpeak = spike_detect_threshold(V0, V1,
		time, tpeak, Vpeak, netpar, runpar, fl)
            else
                println("spike_detect_criterion=",
		runpar.spike_detect_criterion, " should be p or t")
                exit(0)
   	    end

            if spike_detected
	        after_max_vol[jon] = 1
		jpop = netpar.whichpop[jon]
		jon_pop = jon+1-netpar.pinds[jpop]
		
#		push!(rst_ar, [tpeak, jon, jpop, jon_pop])
#               if tpeak > runpar.Tall - runpar.tstat - runpar.epsilon
		    push!(rst_ar[jpop], NamedTuple{(:tspk, :ion_pop)}((tpeak,
		      jon_pop)))
#		end
                runpar.sm ? @printf(fl.rst, "%14.8e %d %d %d %lf\n", tpeak,
		  jon, jpop, jon_pop, Vpeak) : nothing
		if time <= runpar.time_for_recognizing_spikes[jpop]
		    actions_spike_detect(jpop, jon, tpeak, Vpeak, it, time,
                    synstr, spkstr, netpar, runpar, fl)
		end
            end
	else
	     if V0 < runpar.after_min_vol
                 after_max_vol[jon] = 0
             end
	end
    end
    
end #spike_detect

#This function finds the spikes fired by pre-synaptic neurons and updates
#the synaptic variables of the post-synaptic neurons.
function spike_generate(Varbar, Varold, after_max_vol, it, time, rst_ar,
         synstr, spkstr, netpar, runpar, fl)

    spkstr.nsp_in_dt = 0
    spike_detected = false
    jon = 0
    jpop = 0
    tpeak = 0.0
    Vpeak = 0.0

    for kpop in 1:netpar.Npop
        tpop = 500 + (kpop-1) * 1000.0
	if abs(time - tpop) < 1.0e-8
	    jpop = kpop
	    jon = netpar.pinds[jpop]
	    spike_detected = true;
	    tpeak = time;
	    Vpeak = 50.0;
	    printfmtln("generate spikes time={1:f} tpeak={2:f} kpop={3:d}" *
	    " jon={4:d}", time, tpeak, kpop, jon)
	end
    end
	
    if spike_detected
        jon_pop = jon+1-netpar.pinds[jpop]
#	push!(rst_ar, [tpeak, jon, jpop, jon_pop])
  	push!(rst_ar[jpop], NamedTuple{(:tspk, :ion_pop)}((tpeak, jon_pop)))
        runpar.sm ? @printf(fl.rst, "%14.8e %d %d %d %lf\n", tpeak, jon, jpop,
	  jon_pop, Vpeak) : nothing
	actions_spike_detect(jpop, jon, tpeak, Vpeak, it, time, synstr,
          spkstr, netpar, runpar, fl)
    end
    
end #spike_generate

#This function detects a peak of the membrane potential.
function spike_detect_peak(V0, V1, V2, time, tpeak, Vpeak, netpar, runpar, fl)

    if V1 >= V0 && V1 >= V2 && V1 > runpar.Volt_thresh  #detect a spike
        spike_detected = true
        xb = V2 - V0;
        xc = V0 - 2.0 * V1 + V2
        if (fabs(xc) < runpar.epsilon)
            tpeak = time - runpar.deltat
            Vpeak = V1
        else        
            tpeak = time - runpar.deltat + 0.5 * (xb / xc) * runpar.deltat
            Vpeak = V1 - 0.125 * xb * xb / xc
        end
    else
        spike_detected = false, tpeak, Vpeak
    end

  return spike_detected
end #spike_detect_peak

#This function detects a threshold crossing of the membrane potential.
function spike_detect_threshold(V0, V1, time, tpeak, Vpeak, netpar, runpar, fl)

    if V0 >= runpar.Volt_thresh && V1 < runpar.Volt_thresh  #detect a spike
        spike_detected = true
        tpeak = lininter(V1, V0, runpar.Volt_thresh, time - runpar.deltat, time)
        Vpeak = runpar.Volt_thresh
    else
        spike_detected = false
    end

    return spike_detected, tpeak, Vpeak
end #spike_detect_threshold

#This functions makes the necessary actions when a spike is detected.
function actions_spike_detect(jpop, jon, tpeak, Vpeak, it, time, synstr,
    spkstr, netpar, runpar, fl)
    
#   This command stores newly-detected spikes in the array spkstr.sp_in_dt
#   It is used for cases in which the values of synaptic delays tau_delay vary 
#   among synapses.
    spkstr.nsp_in_dt += 1
    spkstr.sp_in_dt[spkstr.nsp_in_dt] = (tspike=tpeak, jpop=jpop, jon=jon)

#   multiple_td_storing_spikes(jpop, jon, tpeak, it, time, synstr, spkstr,
#   netpar, runpar, fl)
    
end #actions_spike_detect

#function multiple_td_storing_spikes(jpop, jon, tpeak, it, time, synstr, spkstr,
#        netpar, runpar, fl)
	 
#   spkstr.nsp_in_dt += 1
##  sp_tuple = NamedTuple{(:tspike, :jpop, :jon)}((tpeak, jpop, jon))
##  spkstr.sp_in_dt[spkstr.nsp_in_dt] = sp_tuple
#   spkstr.sp_in_dt[spkstr.nsp_in_dt] = (tspike=tpeak, jpop=jpop, jon=jon)

#end #multiple_td_storing_spikes

#This function takes all the spikes fired within the time step deltat,      
#adds the specific synaptic delay for each pair of pre- and post-synaptic   
#neurons, and stores the existence of the coming spike to the post-synaptic 
#neuron, at the time that is the sum of spike time + tau_delay, in the      
#array time_spike_delay.                                                    
function multiple_store_spikes_plus_td(it, time, synstr, spkstr, netpar, runpar,         fl)

    for ipop in 1:netpar.Npop, jpop in 1:netpar.Npop
        if netpar.connection_exists[ipop, jpop]
            spkstr.iptr_delay[ipop, jpop] += 1
	    while spkstr.iptr_delay[ipop, jpop] > spkstr.ndelay[ipop, jpop]
	        spkstr.iptr_delay[ipop, jpop] -= spkstr.ndelay[ipop, jpop]
	    end
        end
    end

    for ispk in 1:spkstr.nsp_in_dt
        (tspk, jpop, jon) = spkstr.sp_in_dt[ispk]
	for jcoup in 1:length(netpar.wcoup[jon])
	    ion = netpar.wcoup[jon][jcoup]
	    ipop = netpar.whichpop[ion]
	    tspike_arrive = tspk + netpar.tdcoup[jon][jcoup]
	    idelay_a = trunc(Int64, (runpar.deltat + tspike_arrive - time +
	    runpar.epsilon) / runpar.deltat)
	    idelay_b = idelay_a + spkstr.iptr_delay[ipop, jpop]
	    idelay = idelay_b <= spkstr.ndelay[ipop, jpop] ? 
	      idelay_b : idelay_b - spkstr.ndelay[ipop, jpop]
	    spkstr.time_spike_delay[ion][jpop][idelay] +=
	      netpar.S[ipop, jpop].UU
        end
    end
    
end #multiple_store_spikes_plus_td

# This function stores the number of spikes fired during preBot activity.
function spikes_during_prebot_activity(it, time, synstr, spkstr, netpar,
  runpar, fl)
      
    if synstr.s_ext[1] > runpar.epsilon &&
      synstr.s_ext_previous[1] < runpar.epsilon
        spkstr.prebot_on_for_counting_spikes = true
        spkstr.spk_during_prebot_now = zeros(Int64, netpar.Ntot)
    end

    if spkstr.prebot_on_for_counting_spikes
        for ispk in 1:spkstr.nsp_in_dt
            (tspk, jpop, jon) = spkstr.sp_in_dt[ispk]
	    spkstr.spk_during_prebot_now[jon] += 1
	end
    end

    if synstr.s_ext[1] < runpar.epsilon &&
      synstr.s_ext_previous[1] > runpar.epsilon
        spkstr.prebot_on_for_counting_spikes = false
        spkstr.n_prebot_on_for_counting_spikes += 1
        for ion in 1:netpar.Ntot
	   spkstr.spk_during_prebot[ion] += spkstr.spk_during_prebot_now[ion]
	end
    end
    
end #spikes_during_prebot_activity

#This function compute the decay of synaptic variables.
function decay_post_synaptic_variables(synstr, time, it, netpar, runpar, fl)

    for ion in 1:netpar.Ntot
        ipop = netpar.whichpop[ion]
	for isynvar in 1:synstr.nsynvar[ipop]
	    synstr.synvar[ion][isynvar] *= exp(-runpar.deltat /
	      synstr.synpar_during_run[ipop][isynvar].tsyn);
	end
    end

end #decay_post_synaptic_variables

# This function updates the synaptic variables of post-synaptic neurons in
# response to firing of pre-synaptic neurons. Time delays are considered.
function update_post_synaptic_variables_for_pre_synaptic_spikes(synstr, spkstr,
         time, it, netpar, runpar, fl)

    for ion in 1:netpar.Ntot
        ipop = netpar.whichpop[ion]
	
        for jpop in 1:netpar.Npop
	    if netpar.connection_exists[ipop, jpop]
	        idelay = spkstr.iptr_delay[ipop, jpop]
		if spkstr.time_spike_delay[ion][jpop][idelay] >
		  runpar.epsilon
		    isyn_to_send_1 = synstr.isyn_to_send[ipop][jpop][1]
		    isyn_to_send_2 = synstr.isyn_to_send[ipop][jpop][2]
		    for isynvar in isyn_to_send_1:isyn_to_send_2
		        synstr.synvar[ion][isynvar] += 
	                  spkstr.time_spike_delay[ion][jpop][idelay]
		    end
		    spkstr.time_spike_delay[ion][jpop][idelay] = 0.0    
		end
	    end
	end
    end
    
#   printfmtln("{1:f} {2:f}", time, synstr.synvar[1][2])
end #update_post_synaptic_variables_for_pre_synaptic_spikes

# This function updates the arrays of Vav every time step.
# The arrays are needed to compute the synchrony measure chi.
function update_Vav_arrays(Varbar, Vav, it, time, deltat, netpar, runpar, fl)

    for ipop in 1:netpar.Npop
        Vav[ipop].Vpop = mean(Varbar[netpar.pinds[ipop]:netpar.pinds[ipop+1]-
	  1])[1]
	for ion in netpar.pinds[ipop]:netpar.pinds[ipop+1]-1
	    ion_Vav = ion - netpar.pinds[ipop] + 1
            Vav[ipop].V_avt[ion_Vav]    += Varbar[ion][1] * deltat
            Vav[ipop].V_sq_avt[ion_Vav] += Varbar[ion][1]^2.0 * deltat
	    if Varbar[ion][1] < runpar.spk_threshold
	        Vav[ipop].t_spk_clamp[ion_Vav] += deltat
                Vav[ipop].V_avt_spk_clamp[ion_Vav] += Varbar[ion][1] * deltat
            end
	end
	
        Vav[ipop].Vpop_avt    += Vav[ipop].Vpop   * deltat;
        Vav[ipop].Vpop_sq_avt += Vav[ipop].Vpop^2 * deltat;
    end

end #update_Vav_arrays

# This function updates the arrays of ISIav every time step.
# The arrays are needed to compute CV and CVt.
function update_ISIav_arrays(spkstr, it, time, netpar, runpar, fl)

for ispk in 1:spkstr.nsp_in_dt
    (tspk, jpop, jon) = spkstr.sp_in_dt[ispk]
    spkstr.ISIav.nspk[jon] += 1
    
#   if jon == 1
#       println(" tspk=", tspk, " n=", spkstr.ISIav.nspk[jon])
#   end
    
    if spkstr.ISIav.nspk[jon] >= 2
        delta_tspk = tspk - spkstr.ISIav.tspk_old[jon]
	spkstr.ISIav.av[jon] += delta_tspk
	spkstr.ISIav.avt[jon] += delta_tspk * delta_tspk
	
#       if jon == 1
#           println("dtspk=", delta_tspk, " av=", spkstr.ISIav.av[jon], " avt=",
#  	      spkstr.ISIav.avt[jon])
#       end

        if spkstr.ISIav.nspk[jon] >= 3
	    CVt_now = 2.0* abs(delta_tspk - spkstr.ISIav.delta_tspk_old[jon]) /
	      (delta_tspk + spkstr.ISIav.delta_tspk_old[jon])
	    spkstr.ISIav.CVt[jon] += CVt_now
#           if jon == 1
#               println("CVt_now=", CVt_now, " CVt=", spkstr.ISIav.CVt[jon])
#           end
	end
	
	spkstr.ISIav.delta_tspk_old[jon] =  delta_tspk
    end
    
    spkstr.ISIav.tspk_old[jon] = tspk
    
end

end #update_ISIav_array

#######
#This function computes statistics of spikes, voltages and synaptic values.
function compute_statistics(Varbar, synstr, statstr, netpar, runpar, fl)

    for ipop = 1:netpar.Npop
        statstr.Vav_pop[ipop] = 0.0
        for ion in netpar.pinds[ipop]:netpar.pinds[ipop+1]-1
            statstr.Vav_pop[ipop] += Varbar[ion][1]
	end
	statstr.Vav_pop[ipop] /= netpar.pinds[ipop+1] - netpar.pinds[ipop]
    end

    if runpar.consider_s == 'r'
        for ipop = 1:netpar.Npop
	    statstr.sav_pop[ipop][1] = 0.0
            for ion in netpar.pinds[ipop]:netpar.pinds[ipop+1]-1
                statstr.sav_pop[ipop][1] += Varbar[ion][end]
#	        mean(Varbar[netpar.pinds[ipop]:(netpar.pinds[ipop+1]-1)][4])
	    end
	    statstr.sav_pop[ipop][1] /= netpar.pinds[ipop+1] -
	      netpar.pinds[ipop]
	end
    elseif runpar.consider_s == 's'
        for ipop = 1:netpar.Npop
	    for isynvar in 1:synstr.nsynvar[ipop]
 	        statstr.sav_pop[ipop][isynvar] = 0.0
	        for ion in netpar.pinds[ipop]:netpar.pinds[ipop+1]-1
		    statstr.sav_pop[ipop][isynvar] +=
		      synstr.synvar[ion][isynvar]
		end
		statstr.sav_pop[ipop][isynvar] /= netpar.pinds[ipop+1] -
	          netpar.pinds[ipop]
	    end
        end
    end

end #compute_statistics

#This function computes statistics of synaptic values.
function compute_syn_statistics(statstr, netpar, runpar, fl)

    Sp = statstr.sav_pop[1][2] + statstr.sav_pop[2][1]
    Sm = statstr.sav_pop[1][2] - statstr.sav_pop[2][1]
    statstr.Sp_av  += Sp * runpar.deltat
    statstr.Sp_avt += Sp * Sp * runpar.deltat
    statstr.Sm_av  += Sm * runpar.deltat
    statstr.Sm_avt += Sm * Sm * runpar.deltat
#   println(Sp, " ", Sm, " ", statstr.Sp_av, " ", statstr.Sp_avt)
    
end #compute_syn_statistics

#This function writes the data as functions of time
function pr_fct(Varbar, statstr, synstr, netpar, runpar, time, it, fl)

    icol_wrt = 1
    
#   write(fl.col, "$time")
    runpar.first_print ? printfmtln(fl.out, "\ncol numbers\ntime {1:d}",
      icol_wrt) : nothing
    @printf(fl.col, " %.5lf", time)

    for ipop = 1:netpar.Npop
	icol_wrt += 1
	runpar.first_print ? printfmtln(fl.out, "{1:d} Vav_pop[{2:d}]",
	  icol_wrt, ipop) : nothing
        printfmt(fl.col, " {1:f}", statstr.Vav_pop[ipop])
	
	if runpar.consider_s == 'r'
	    icol_wrt += 1
	    runpar.first_print ? printfmtln(fl.out, "{1:d} sav_pop", icol_wrt) :
	      nothing
            printfmt(fl.col, " {1:f}", statstr.sav_pop[ipop][1])
#           printfmt(fl.col, " {1:f}", synstr.savr_pop[ipop])
         elseif runpar.consider_s == 's'
	    for isynvar in 1:synstr.nsynvar[ipop]
	        icol_wrt += 1
	        if runpar.first_print 
	             printfmtln(fl.out, "{1:d} sav_pop[{2:d}][{3:d}]",
		     icol_wrt, ipop, isynvar)
		end
                printfmt(fl.col, " {1:f}", statstr.sav_pop[ipop][isynvar])
#		icol_wrt += 1
#		if runpar.first_print		
#		    [printfmtln(fl.out, " {1:d}", synstr.synvar[ion][1]) for ion in 1:2]
#                end
#                for iwrite in 1:2
#	            ion = netpar.pinds[ipop] - 1 + iwrite
#		    printfmt(fl.col, " {1:f}", synstr.synvar[ion][1])
#		end
### to detete
#               ion = netpar.pinds[ipop] # first neuron in a population
#		icol_wrt += 1
#		runpar.first_print ?
#		printfmtln(fl.out, "{1:d} synvar[{2:d}][{3:d}]", icol_wrt, ion,
#		  isynvar) : nothing
#	        printfmt(fl.col, " {1:f}", synstr.synvar[ion][isynvar])
###		
	    end
	else
	    println("consider_s =", runpar.consider_s, " should be r or s!")
	    exit(0)
        end

        ion = netpar.pinds[ipop] # first neuron in a population
	icol_wrt += 1
	runpar.first_print ? printfmtln(fl.out, "{1:d} Isyn_cont[{2:d}]",
	  icol_wrt, ion) : nothing
	printfmt(fl.col, " {1:f}", synstr.Isyn_cont[ion] )

        for iwrite = 1:length(runpar.nwrite[ipop])
	    ion = netpar.pinds[ipop] - 1 + runpar.nwrite[ipop][iwrite]
   	    icol_wrt += 1
	    runpar.first_print ? printfmtln(fl.out, "{1:d} Varbar[{2:d}]",
	      icol_wrt, ion) : nothing
            printfmt(fl.col, " {1:f}", Varbar[ion][1])
 	end
    end
    write(fl.col, "\n")


#   ion_trj = 1
    ipop_trj = netpar.whichpop[runpar.ion_trj]
#   ipop_trj = 1
#   ion_trj = netpar.pinds[ipop_trj]
    printfmt(fl.trj, "{1:10.5f} {2:12.8f}", time, Varbar[runpar.ion_trj][1])
    for ieq in 2:netpar.C[ipop_trj].Neq
        printfmt(fl.trj, " {1:10.8f}", Varbar[runpar.ion_trj][ieq])
    end
    printfmt(fl.trj, "\n")
#   println(fl.trj, time, " ", Varbar[runpar.ion_trj][1], " ",
#   Varbar[runpar.ion_trj][2], " ", Varbar[runpar.ion_trj][3], " ",
#   Varbar[runpar.ion_trj][4])

    runpar.first_print ? printfmtln(fl.out, " ") : nothing
    runpar.first_print = false
end #pr_fct


# This function compute voltage statistics based on the data obtained during
# the process of solving the differential euqations.
function compute_voltage_statistics(Vav, statstr, netpar, runpar, av, fl)

    for ipop in 1:netpar.Npop
        sigma_V_sq = zeros(netpar.Ncells[ipop])
	
        Vav[ipop].Vpop_avt    /= runpar.tstat
        Vav[ipop].Vpop_sq_avt /= runpar.tstat
        sigma_Vpop_sq = Vav[ipop].Vpop_sq_avt - Vav[ipop].Vpop_avt^2
	
	pop_av_sigma_V_sq = 0.0
	non_with_V_below_threshold = 0
	
	for ion in 1:netpar.Ncells[ipop]
	     Vav[ipop].V_avt[ion]    /= runpar.tstat
             Vav[ipop].V_sq_avt[ion] /= runpar.tstat
             sigma_V_sq[ion] = Vav[ipop].V_sq_avt[ion] -
	       Vav[ipop].V_avt[ion]^2.0
             pop_av_sigma_V_sq += sigma_V_sq[ion]

             if Vav[ipop].t_spk_clamp[ion] > runpar.epsilon
	         non_with_V_below_threshold += 1
	         Vav[ipop].V_avt_spk_clamp[ion] /= Vav[ipop].t_spk_clamp[ion]
		 Vav[ipop].Vpop_spk_clamp += Vav[ipop].V_avt_spk_clamp[ion]
             else
	         Vav[ipop].V_avt_spk_clamp[ion] = -999.9
	     end
	end
	pop_av_sigma_V_sq /= netpar.Ncells[ipop]

        if pop_av_sigma_V_sq <= 0.0
            Vav[ipop].chi = -9999.0
        elseif sigma_Vpop_sq < 0.0
            Vav[ipop].chi = -9998.0
        else
            Vav[ipop].chi = sqrt(sigma_Vpop_sq / pop_av_sigma_V_sq)
        end

        if non_with_V_below_threshold > 0
	    Vav[ipop].Vpop_spk_clamp /= non_with_V_below_threshold
	else
	    Vav[ipop].Vpop_spk_clamp = -997.7
	end

        printfmtln(fl.out, "ipop={1:f}", ipop)
	printfmtln(fl.out,
	  "Vpop_avt={1:f} Vpop_sq_avr={2:f} sigma_Vpop_sq={3:f}",
          Vav[ipop].Vpop_avt, Vav[ipop].Vpop_sq_avt, sigma_Vpop_sq)
        printfmtln(fl.out, "pop_av_sigma_V_sq={1:f} chi={2:f}",
	  pop_av_sigma_V_sq, Vav[ipop].chi)
 	  
	av[ipop].V_avt = Vav[ipop].Vpop_avt
	av[ipop].chi = Vav[ipop].chi
	av[ipop].V_avt_spk_clamp = Vav[ipop].Vpop_spk_clamp
    end
    
    statstr.Sp_av  /= runpar.tstat
    statstr.Sp_avt /= runpar.tstat
    statstr.Sm_av  /= runpar.tstat
    statstr.Sm_avt /= runpar.tstat

    statstr.Sp_sd = sd_cal(statstr.Sp_av, statstr.Sp_avt, runpar, fl)
    statstr.Sm_sd = sd_cal(statstr.Sm_av, statstr.Sm_avt, runpar, fl)
    printfmtln(fl.out, "Sp_av={1:f} Sp_sd={2:f} Sm_av={3:f} Sm_sd={4:f}",
      statstr.Sp_av, statstr.Sp_sd,statstr.Sm_av, statstr.Sm_sd)
   
    av[1].Sp_av = statstr.Sp_av
    av[1].Sp_sd = statstr.Sp_sd
    av[1].Sm_av = statstr.Sm_av
    av[1].Sm_sd = statstr.Sm_sd

end #compute_voltage_statistics

# This function compute spike statistics based on the data obtained during
# the process of solving the differential euqations.
function compute_spike_statistics(ISIav::ISI_aver, spkstr, netpar, runpar, av,
  fl)

    for ion in 1:netpar.Ntot
        if ISIav.nspk[ion] <= 2
	    ISIav.CV[ion] = -999.9
	    ISIav.CVt[ion] = -999.8
	else
	    ISIav.av[ion] /= ISIav.nspk[ion]-1
	    ISIav.avt[ion] /= ISIav.nspk[ion]-1
	    
	    diff = ISIav.avt[ion] - ISIav.av[ion] * ISIav.av[ion]
	    if diff >= 0.0
	        sd = sqrt(diff)
		ISIav.CV[ion] = sd / ISIav.av[ion]
            elseif diff >= -runpar.epsilon
	        sd = 0.0
		ISIav.CV[ion] = 0.0
	    else
	        println("diff=", diff, " <= 0!")
		exit(0)
	    end
	    
	    if ISIav.nspk[ion] <= 3
	        ISIav.CVt[ion] = -999.85
	    else
	        ISIav.CVt[ion] /= ISIav.nspk[ion]-2
	    end
	end
    end

    println(fl.out, "  ion    CV      CVt      av")
    for ion in 1:netpar.Ntot
        printfmtln(fl.out, "{1:d} {2:f} {3:f} {4:f}", ion, ISIav.CV[ion],
	  ISIav.CVt[ion], ISIav.av[ion])
    end

#   Compute population-averaged values
    for ipop in 1:netpar.Npop
        ISIav.n_silent[ipop] = 0
	
	ISIav.CV_pop_av[ipop]   = 0.0
	ISIav.CVt_pop_av[ipop] = 0.0
	ISIav.T_pop_av[ipop]   = 0.0
	
        for ion in netpar.pinds[ipop]:netpar.pinds[ipop+1]-1
	    if ISIav.nspk[ion] <= 3
	        ISIav.n_silent[ipop] += 1
	    else
	        ISIav.CV_pop_av[ipop] += ISIav.CV[ion]
	        ISIav.CVt_pop_av[ipop] += ISIav.CVt[ion]
	        ISIav.T_pop_av[ipop] += ISIav.av[ion]
	    end
	end

        if ISIav.n_silent[ipop] == netpar.Ncells[ipop]
#           all silent
            ISIav.CV_pop_av[ipop]  = -998.0
	    ISIav.CVt_pop_av[ipop] = -998.0
	    ISIav.T_pop_av[ipop]   = -998.0
	else
	     ISIav.CV_pop_av[ipop]  /= netpar.Ncells[ipop] -
	       ISIav.n_silent[ipop]
	     ISIav.CVt_pop_av[ipop] /= netpar.Ncells[ipop] -
	       ISIav.n_silent[ipop]
	     ISIav.T_pop_av[ipop]   /= netpar.Ncells[ipop] -
	       ISIav.n_silent[ipop]
	end
	
	av[ipop].n_silent   = ISIav.n_silent[ipop]
	av[ipop].CV_pop_av  = ISIav.CV_pop_av[ipop]
	av[ipop].CVt_pop_av = ISIav.CVt_pop_av[ipop]
 	av[ipop].T_pop_av   = ISIav.T_pop_av[ipop]
    end
    
    println(fl.out, "ipop  n_silent  CV  CVt  T_av")
    for ipop in 1:netpar.Npop
        printfmtln(fl.out, "{1:d} {2:d} {3:f} {4:f} {5:f}", ipop,
	  ISIav.n_silent[ipop], ISIav.CV_pop_av[ipop], ISIav.CVt_pop_av[ipop],
	    ISIav.T_pop_av[ipop])
    end

    println(fl.out, "\nspk_during_prebot")
    spk_during_prebot_avt_avpop = zeros(netpar.Npop)
    ncells_spk_during_prebot_avt_avpop = zeros(Int64, netpar.Npop)
    spk_during_prebot_avt = zeros(netpar.Ntot)
        
    for ipop = 1:netpar.Npop
	non_min = netpar.pinds[ipop]
	non_max = netpar.pinds[ipop+1]-1
        for ion in non_min:non_max
	    if spkstr.n_prebot_on_for_counting_spikes > 0
                spk_during_prebot_avt[ion] = 1.0 *
	          spkstr.spk_during_prebot[ion] /
	          (spkstr.n_prebot_on_for_counting_spikes)
		ncells_spk_during_prebot_avt_avpop[ipop] += 1
	        spk_during_prebot_avt_avpop[ipop] += spk_during_prebot_avt[ion]
            end
	end
	
	if ncells_spk_during_prebot_avt_avpop[ipop] > 0
   	    spk_during_prebot_avt_avpop[ipop] /= netpar.Ncells[ipop]
	else
	    spk_during_prebot_avt_avpop[ipop] = -999.9
	end
	
	printfmtln(fl.out, "ipop={1:d} spk_during_prebot_avt_avpop={2:f}",
	  ipop, spk_during_prebot_avt_avpop[ipop])
    end

    for ipop = 1:netpar.Npop
        av[ipop].spk_during_prebot_avt_avpop = spk_during_prebot_avt_avpop[ipop]
    end

end #compute_spike_statistics

# This function Writes parameters on parthohop, to be transferred to hop_cal.jl.
function fill_in_par_to_hop(sval, netpar, runpar, partohop, fl)
    partohop.scan_type = sval.scan_type
    partohop.open_file_for_par = sval.open_file_for_par
    partohop.Tall = runpar.Tall
    partohop.determine_time_prebot = runpar.determine_time_prebot
    partohop.Tper_syn_ext = netpar.Tper_syn_ext
    partohop.Tup_syn_ext = netpar.Tup_syn_ext
    partohop.ipop_cal_wsk = runpar.ipop_cal_wsk
    partohop.angle_cal = runpar.angle_cal
    partohop.whisk_cal = runpar.whisk_cal
    partohop.nonF = netpar.Ncells[runpar.ipop_cal_wsk]
    partohop.tstat = runpar.tstat
    partohop.frac_sd_mnmx = runpar.frac_sd_mnmx
    T_diff_pre_interval = deepcopy(runpar.T_diff_pre_interval)
    partohop.whisker = deepcopy(netpar.whisker)
end


# This function calculates the standard deviation
function sd_cal(xav, xtav, runpar, fl)

    sd = 0.0
    diff = xtav - xav * xav
    if diff >= 0.0
        sd = sqrt(diff)
    elseif diff > -runpar.epsilon
        sd = 0.0
    else
        sd = -999.9
    end

    return sd
end #sd_cal

#Linear interpolation 
function lininter(x1,  x2, xc, y1, y2)

    linter = ((xc-x1)*y2+(x2-xc)*y1) / (x2-x1)
    return linter
end  #lininter
