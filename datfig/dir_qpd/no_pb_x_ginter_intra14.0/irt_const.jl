sval = scan_val()
sval.scan_type = 'y'
sval.parname = ["netpar.S_orig[1,2].syn_receptor_par_ar[1].gsyn",
"netpar.S_orig[2,1].syn_receptor_par_ar[1].gsyn = netpar.S_orig[1,2].syn_receptor_par_ar[1].gsyn"]
sval.parname_gr = "g\\ssyn\\N\\S\\h{-1.5}I1I1\\N"
sval.parmin = 20.0
sval.parmax = 46.0
sval.npar = 26
sval.nrepeat = 5
sval.seed = 9751
sval.open_file_for_par = false

svala = scan_val()
svala.parname = ["syncoupparII.gel"]
svala.parname_gr =  "g\\sel\\N\\S\\h{-0.8}II\\N"
svala.parmin = 0.02
svala.parmax = 1.02
svala.npar = 20

#FN
cellparF = cell_par()
cellparF.Neq = 7   # (V, h, n, z, u, r, ss)
cellparF.func_name = "F"

cellparF.Cms = 1.0
cellparF.gNa = 100.0 
cellparF.gNap = 0.04
cellparF.gKdr = 20.0
cellparF.gKm  = 1.0 * 0.3
cellparF.del_gKm = 0.0
cellparF.gAHP = 10.0 * 0.0
cellparF.gh   = 0.05
cellparF.gL   = 0.12
cellparF.del_gL = 0.06

cellparF.VNa  =  55.0
cellparF.VK   = -90.0
cellparF.Vh   = -27.4
cellparF.VL   = -70

cellparF.thetam = -28.0
cellparF.sigmam = -7.8
cellparF.thetah = -50.0
cellparF.sigmah =  7.0
cellparF.thetan = -23.0
cellparF.sigman = -15.0
cellparF.thetap = -53.0
cellparF.sigmap = -5.0
cellparF.thetaz = -45.0
cellparF.sigmaz = -4.25
cellparF.tauz   = 75.0
cellparF.thetau = -25.0
cellparF.sigmau = -3
cellparF.tauu   = 75.0
cellparF.thetar = -83.9
cellparF.sigmar =  7.4

cellparF.Iapp = 3.1
cellparF.phi  = 1.0

cellparF.gexc_const = 0.0 # 0.05
cellparF.Vexc = 0.0

cellparF.ginh_const = 0.0

cellparF.gsyn_ext = 0.0
cellparF.del_gsyn_ext = 0.0
cellparF.Vsyn_ext = -80.0

cellparF.Vinc1 = -65.0
cellparF.Vinc2 = -55.0

cellparF.alpha = 12.0
cellparF.beta = 0.1
cellparF.theta_syn = 0.0
cellparF.sigma_syn = 2.0

cellparF.initial_val = [-35.0, 0.4, 0.1, 0.1, 0.1, 0.1, 0.1]

#IRT - based on F (motorneuron, mn) with two slow variables
cellparI_mn = deepcopy(cellparF)
cellparI_mn.Neq = 6   # (V, h, n, z, s, ss)
cellparI_mn.func_name = "I_mn"

cellparI_mn.gNa = 100.0
cellparI_mn.gNap = 0.04
cellparI_mn.gKdr = 20.0
cellparI_mn.gKm = 7.0
cellparI_mn.del_gKm = 3.0
cellparI_mn.gCas = 0.0
cellparI_mn.gh   = 0.0
cellparI_mn.gL   = 0.12
cellparI_mn.del_gL = 0.06

cellparI_mn.VCa = 120

cellparI_mn.thetaz = -28.0
cellparI_mn.sigmaz = -3.0
cellparI_mn.tauz   = 83.0

cellparI_mn.thetas = 12.0
cellparI_mn.sigmas = -12.0
cellparI_mn.taus = 4.5

cellparI_mn.Iapp = 20.0

cellparI_mn.gexc_const = 0.0
cellparI_mn.Vexc = 0.0

cellparI_mn.ginh_const = 0.0

cellparI_mn.gsyn_ext = 0.0
cellparI_mn.del_gsyn_ext = 0.0
cellparI_mn.Vsyn_ext = -80.0

cellparI_mn.Vinc1 = -65.0
cellparI_mn.Vinc2 = -55.0

cellparI_mn.initial_val = [-52.0, 0.4, 0.1, 0.1, 0.1, 0.1]

cellparI1 = deepcopy(cellparI_mn)

cellparI2 = deepcopy(cellparI_mn)
cellparI2.gsyn_ext = 0.0
cellparI2.del_gsyn_ext = 0.0

#-------------------
#IRT1 <- IRT1
syncoupparI1I1 = syn_coup_par()
syncoupparI1I1.syn_receptor_par_ar = Array{syn_receptor_par,1}(undef,1)

#GABAA
syncoupparI1I1.syn_receptor_par_ar[1] = syn_receptor_par()
syncoupparI1I1.syn_receptor_par_ar[1].gsyn = 14.0
syncoupparI1I1.syn_receptor_par_ar[1].Vsyn = -80.0
syncoupparI1I1.syn_receptor_par_ar[1].tsynr = 0.0
syncoupparI1I1.syn_receptor_par_ar[1].tsynd = 10.0
syncoupparI1I1.syn_receptor_par_ar[1].is_it_NMDA = false

syncoupparI1I1.Kin = 25.0
syncoupparI1I1.UU = 1.0
syncoupparI1I1.tau_delay = 0.2 * 0.0
syncoupparI1I1.Del_tau_delay = 0.1 * 0.0

syncoupparI1I1.gel = 0.0
syncoupparI1I1.Kel = 25.0

#IRT1 <- IRT2
syncoupparI1I2 = deepcopy(syncoupparI1I1)
syncoupparI1I2.syn_receptor_par_ar[1].gsyn = 6.0
#syncoupparI1I2.gel = 0.0

#IRT2 <- IRT1
syncoupparI2I1 = deepcopy(syncoupparI1I2)
syncoupparI2I1.syn_receptor_par_ar[1].gsyn = 6.0
#syncoupparI2I1.gel = 0.0

#IRT2 <- IRT2
syncoupparI2I2 = deepcopy(syncoupparI1I1)
syncoupparI2I2.syn_receptor_par_ar[1].gsyn = 14.0
#syncoupparI2I2.gel = 0.0

#FN <- IRT1
syncoupparFI1 = deepcopy(syncoupparI1I1)

syncoupparFI1.syn_receptor_par_ar[1].gsyn = 3.0
syncoupparFI1.syn_receptor_par_ar[1].Vsyn = -80.0
syncoupparFI1.Kin = 25.0

syncoupparFI1.tau_delay = 0.6 * 0.0
syncoupparFI1.Del_tau_delay = 0.1 * 0.0

syncoupparFI1.gel = 0.0

#FN <- IRT2
syncoupparFI2 = deepcopy(syncoupparFI1)
syncoupparFI2.syn_receptor_par_ar[1].gsyn = 0.0
syncoupparFI2.Kin = 25.0

#Whisker
whiskerstr = whisker_par()
whiskerstr.r0 = 1.9
whiskerstr.taur = 5.0
whiskerstr.tauc = 6.0
whiskerstr.tauw = 20.0
whiskerstr.AA = 100.0

#Network
netpar = net_par()
netpar.Ncells_orig = [100, 100, 100] #number of cells in each population: [IRT]
netpar.Npop = length(netpar.Ncells_orig)
netpar.inhibition_factor = 1.0
netpar.inhibition_factor_consider = ["I1"]
netpar.Kfactor = 0.0

netpar.ar_cell_str = ["I1", "I2", "F"]
netpar.C_orig = [cellparI1, cellparI2, cellparF]
netpar.C = deepcopy(netpar.C_orig)

#                            I1   I1    F 
netpar.connection_exists = [true true false;   # I1
                            true true false;   # I2
			    true true false]   # F

netpar.S_orig = Array{syn_coup_par}(undef, 3, 3)
netpar.S_orig[1, 1] = syncoupparI1I1
netpar.S_orig[1, 2] = syncoupparI1I2

netpar.S_orig[2, 1] = syncoupparI2I1
netpar.S_orig[2, 2] = syncoupparI2I2

netpar.S_orig[3, 1] = syncoupparFI1
netpar.S_orig[3, 2] = syncoupparFI2

netpar.S = deepcopy(netpar.S_orig)

netpar.Kel_scale = 25.0

netpar.rho_concur = 1.0
netpar.noise = 0.0

netpar.type_syn_ext = 'r' # p - periodic, r - random.
netpar.Tper_syn_ext = 700.0
netpar.Trand_syn_ext = 150.0
netpar.Tup_syn_ext = 70.0

netpar.whisker = whiskerstr

#Parameters of simulation running
runpar = run_par()

runpar.epsilon = 1.0e-11
runpar.seed = sval.seed
runpar.Tall = 6000.0
runpar.deltat = 0.01
runpar.incond = 'b' # r - read initial conditions, b - auxiliary variables
                    # calculated from V, s - steady state, a - for checking.
runpar.Vincond_add = 0.01
runpar.Vincond_rand = 0.01
runpar.NT = round(Int64, runpar.Tall/runpar.deltat)
runpar.method = 'r'
runpar.Volt_thresh = -20.0
runpar.spike_detect_criterion = 't'
runpar.after_min_vol = -33.0
runpar.consider_s = 's' # variable s: r - pre-synaptic, s - post-synaptic.
runpar.generate_one_psp = 'n'  # n - network simulation, y - only single PSPs.

runpar.nwrite=[[1, 8],  #I1
               [1, 2],  #I2
	       [1]]     #F
runpar.twrite = 2
runpar.tmcol = 7000.0
runpar.tstat = 5000.0
runpar.traster = 6000.0
runpar.ion_trj = 1
runpar.smforce = 'l' # p - always print, n - always no print, l - leave as is.
runpar.sp = 1
runpar.time_for_recognizing_spikes = [1.0e12, 1.0e12, 1.0e12]
runpar.determine_time_prebot = 'e' # o - one neuron, e - external gsyn.
runpar.spk_threshold = -40.0
runpar.ipop_cal_wsk = 3 #F
runpar.angle_cal = 'p' # l - linear, p - power
runpar.whisk_cal = true
runpar.T_diff_pre_interval = [40.0, 140.0]
runpar.frac_sd_mnmx = 0.8

runpar.dt_bin = 1.0
runpar.lags = 1500
runpar.lag_for_Tper = 6
