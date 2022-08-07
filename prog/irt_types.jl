mutable struct scan_val
    parmin::Float64
    parmax::Float64
    par::Float64
    seed::Int64
    npar::Int64
    ipar::Int64
    nrepeat::Int64
    irepeat::Int64
    scan_type::Char
    parname::Array{String,1}
    parname_gr::String
    open_file_for_par::Bool
    
    scan_val() = new()
end

mutable struct avr_val
    Npop::Int64
    pop_name::Char
    V_avt::Float64
    chi::Float64
    V_avt_spk_clamp::Float64
    n_silent::Float64
    CV_pop_av::Float64
    CVt_pop_av::Float64
    T_pop_av::Float64

    Sp_av::Float64
    Sp_sd::Float64
    Sm_av::Float64
    Sm_sd::Float64
    
    Tper::Float64
    Mfr::Float64
    mode::Int64
    Tper_from_cc_av::Float64
    CV_bur_pop_av::Float64
    CVt_bur_pop_av::Float64

    spk_during_prebot_avt_avpop::Float64
    
    avr_val() = new()
end

mutable struct cell_par
    Neq::Int64
    func_name::String
    
    Cms::Float64
    gNa::Float64
    gNap::Float64
    gKdr::Float64
    gKs::Float64
    gKm::Float64
    del_gKm::Float64
    gAHP::Float64
    gh::Float64
    gCas::Float64
    gL::Float64
    del_gL::Float64

    VNa::Float64
    VK::Float64
    Vh::Float64
    VCa::Float64
    VL::Float64

    thetam::Float64
    sigmam::Float64
    thetap::Float64
    sigmap::Float64
    thetah::Float64
    sigmah::Float64
    thetan::Float64
    sigman::Float64
    taunb0::Float64
    taunb::Float64
    thetaz::Float64
    sigmaz::Float64
    tauz::Float64
    thetau::Float64
    sigmau::Float64
    tauu::Float64
    thetar::Float64
    sigmar::Float64
    thetaks::Float64
    sigmaks::Float64
    tauksb0::Float64
    tauksb::Float64

    thetas::Float64
    sigmas::Float64
    taus::Float64

    Iapp::Float64
    phi::Float64
 
    gexc_const::Float64
    Vexc::Float64
    ginh_const::Float64
    
    gsyn_ext::Float64
    del_gsyn_ext::Float64
    Vsyn_ext::Float64
    
    Vinc1::Float64
    Vinc2::Float64
    
    alpha::Float64
    beta::Float64
    theta_syn::Float64
    sigma_syn::Float64

    initial_val::Array{Float64,1}

    steady_state_var::Function
    update_cell!::Function
    cell_par() = new()
end

mutable struct syn_receptor_par
    gsyn::Float64
    Vsyn::Float64
    tsynd::Float64
    tsynr::Float64
    thetanp::Float64
    sigmanp::Float64
    is_it_NMDA::Bool
    
    syn_receptor_par() = new()
end

mutable struct syn_coup_par
    receptor_type::Char
    syn_receptor_par_ar::Array{syn_receptor_par,1}
    Kin::Float64
    UU::Float64
    tau_delay::Float64
    Del_tau_delay::Float64
    
    gel::Float64
    Kel::Float64    
    gel_one_psp::Float64


    syn_coup_par() = new()
end

mutable struct whisker_par
    r0::Float64
    taur::Float64
    tauc::Float64
    tauw::Float64
    AA::Float64

    whisker_par() = new()
end

mutable struct net_par
    Ncells_orig::Array{Int64,1}  #number of cells in each population: [IRT, PreBot, FN]
    Ncells::Array{Int64,1}      #number of cells in each population: [IRT, PreBot, FN]
    Ntot::Int64
    Npop::Int64
    pinds::Array{Int64,1}
    whichpop::Array{Int64,1}    
    ar_cell_str::Array{String,1}
    C::Array{cell_par,1}
    C_orig::Array{cell_par,1}
    S::Array{syn_coup_par,2}
    S_orig::Array{syn_coup_par,2}
#   Neq_max::Int64
    connection_exists::Array{Bool,2}
    Kel_scale::Float64
    rho_concur::Float64
    noise::Float64
    inhibition_factor::Float64
    inhibition_factor_consider::Array{String,1}
    Kfactor::Float64
    
    type_syn_ext::Char
    Tper_syn_ext::Float64
    Trand_syn_ext::Float64
    Tup_syn_ext::Float64

    gKmar::Array{Float64,1}
    gLar::Array{Float64,1}
    gsyn_ext_ar::Array{Float64,1}

    nwcoup::Array{Array{NamedTuple{(:ipop, :ion_min, :ion_max),
            Tuple{Int64, Int64, Int64}}}}
    wcoup::Array{Array{Int64,1},1}
    tdcoup::Array{Array{Float64,1}}

#   necoup::Array{Array{Int64,1}
    ecoup::Array{Array{Tuple{Int64, Float64},1},1}   # [ion][jon][2]

    whisker::whisker_par  

    net_par() = new()
end

mutable struct run_par
    epsilon::Float64
    seed::Int64
    Tall::Float64
    deltat::Float64
    NT::Int64
    incond::Char
    Vincond_add::Float64
    Vincond_rand::Float64
    method::Char
    Volt_thresh::Float64
    spike_detect_criterion::Char
    after_min_vol::Float64
    consider_s::Char
    generate_one_psp::Char
    nwrite::Array{Array{Int64,1},1}
    twrite::Int64
    first_print::Bool
    tmcol::Float64
    tstat::Float64
    traster::Float64
    ion_trj::Int64
    sm::Bool
    sp::Bool
    smforce::Char
    time_for_recognizing_spikes::Array{Float64,1}
    determine_time_prebot::Char
    spk_threshold::Float64
    ipop_cal_wsk::Int64
    angle_cal::Char
    whisk_cal::Bool
    open_file_for_par::Bool
    T_diff_pre_interval::Array{Float64,1}
    frac_sd_mnmx::Float64
    
    dt_bin::Float64
    lags::Int64
    lag_for_Tper::Int64
    
    
    run_par() = new()
end

mutable struct syn_par_during_run
    gsyn_one_psp::Float64
    Vsyn::Float64
    tsyn::Float64
    is_it_NMDA::Bool
    
    syn_par_during_run() = new()
end

mutable struct syn_str
    nsynvar::Array{Int64,1}                                 #[Npop]
    synpar_during_run::Array{Array{syn_par_during_run,1},1} #[ipop][isynvar]
    isyn_to_send::Array{Array{Array{Int64,1},1},1}          #[ipop][jpop][2]
    synvar::Array{Array{Float64,1},1}                       #[ion][isynvar]

    savr_pop::Array{Float64,1}
    Isyn_cont::Array{Float64,1}

    s_ext::Array{Int64,1}                                   #[Npop]
    s_ext_previous::Array{Int64,1}
    Tper_syn_ext_rand::Float64
    T_start_cycle::Float64

    Iel_cont::Array{Float64,1}
    syn_str() = new()
end

mutable struct stat_str
    Vav_pop::Array{Float64,1}
    sav_pop::Array{Array{Float64,1}}
    Sp_av::Float64
    Sp_avt::Float64
    Sp_sd::Float64
    Sm_av::Float64
    Sm_avt::Float64
    Sm_sd::Float64
    
    stat_str() = new()
end

mutable struct V_aver
    Vpop::Float64
    chi::Float64
    Vpop_spk_clamp::Float64

    Vpop_avt::Float64
    Vpop_sq_avt::Float64
    V_avt::Array{Float64,1}
    V_sq_avt::Array{Float64,1}
    t_spk_clamp::Array{Float64,1}
    V_avt_spk_clamp::Array{Float64,1}
    V_aver() = new()
end

mutable struct ISI_aver
    nspk::Array{Int64,1}
    tspk_old::Array{Float64,1}
    delta_tspk_old::Array{Float64,1}
    av::Array{Float64,1}
    avt::Array{Float64,1}
    CV::Array{Float64,1}
    CVt::Array{Float64,1}
    
    n_silent::Array{Int64,1}
    CV_pop_av::Array{Float64,1}
    CVt_pop_av::Array{Float64,1}
    T_pop_av::Array{Float64,1}
    ISI_aver() = new()
end

mutable struct spk_str
    nsp_in_dt::Int64
    sp_in_dt::Array{NamedTuple{(:tspike, :jpop, :jon), Tuple{Float64, Int64,
              Int64}}}

    ndelay::Array{Int64,2}
    iptr_delay::Array{Int64,2}
    
    time_spike_delay::Array{Array{Array{Float64,1},1},1} # [ion, jpop, idelay]

    Vav::Array{V_aver,1}
    ISIav::ISI_aver

#	      iptr_delay[ipop][jpop]
#             ndelay[ipop][jpop]
#              >time_spike_delay[ipop][jpop][ion][idelay]

    spk_during_prebot::Array{Int64,1}
    spk_during_prebot_now::Array{Int64,1}
    prebot_on_for_counting_spikes::Bool
    n_prebot_on_for_counting_spikes::Int64
    spk_str() = new()
end

mutable struct rst_pb_str
    rst_ar::Array{Array{NamedTuple{(:tspk, :ion_pop), Tuple{Float64, Int64}},1}}
    prebot_burst_start::Array{Float64,1}

    rst_pb_str() = new()
end

mutable struct bur_str
    Tper::Float64
    Mfr::Array{Float64,1}
    mode::Int64
    Tper_from_cc_av::Float64
    nbin::Int64
    
    bur_str() = new()
end

mutable struct fl_st
    col::IOStream
    trj::IOStream
    sex::IOStream
    wsk::IOStream
    rst::IOStream
    avr::IOStream
    avc::IOStream
    res::IOStream
    hei::IOStream
    dtm::IOStream
    tpp::IOStream
    aws::IOStream
    bur::IOStream
    out::IOStream
    
    fl_st() = new()
end
