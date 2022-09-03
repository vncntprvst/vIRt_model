sval = scan_val()
sval.scan_type = 'y'
sval.parname = ["netpar.I0"]
sval.parname_gr = "g\\ssyn\\N\\S\\h{-1.5}I1I1\\N"
sval.parmin = 0.0
sval.parmax = 30.0
sval.npar = 120
sval.nrepeat = 1
sval.seed = 9751

svala = scan_val()
svala.parname = ["netpar.IrB"]
svala.parname_gr =  "g\\sel\\N\\S\\h{-0.8}II\\N"
svala.parmin = 0.02
svala.parmax = 1.02
svala.npar = 20

#Network
netpar = net_par()
netpar.Npop = 3
netpar.initial_val = [0.1, 0.0, 0.15, 0.3,  0.01, 0.02,  0.0]

netpar.I0_threshold = 0.29
netpar.I0 = 20.0
netpar.taua = 83.0
netpar.taus = 10.0 #0.01 #10.0
netpar.beta = 0.01746  # 0.024
netpar.gamma = 0.02471  #0.023
netpar.gKm = 7.0

netpar.grB = 0.0
netpar.TB = 700.0
netpar.Tup = 70.0

netpar.I0_threshold_F = 0.46
netpar.I0_F = 3.1
netpar.taua_F = 75.0
netpar.beta_F = 0.0305 #0.0222
netpar.gamma_F = 0.0617 #0.05367
netpar.gKm_F  = 0.3

netpar.DelV = 27.0 #y18.0
netpar.DelV_F = 20.0 #y18.0
netpar.gintra = 0.0
netpar.ginter = 6.0

netpar.gFr = 3.0

netpar.inhibition_factor = 1.0

netpar.tauw = 20.0

#Parameters of simulation running
runpar = run_par()

runpar.epsilon = 1.0e-10
runpar.Tall = 16000.0
runpar.deltat = 0.02
runpar.method = 'r'
runpar.incond = 'r' # r - read initial conditions.
runpar.tstat = 14000.0
runpar.twrite = 1
runpar.tmcol = 7000.0
runpar.smforce = 'l' # p - always print, n - always no print, l - leave as is.
runpar.n_hist_store = 7 # odd number
runpar.open_file_for_par = false
runpar.min_theta_for_osci = 1.0
