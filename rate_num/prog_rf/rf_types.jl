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
    Mav::Array{Float64,1}
    Msd::Array{Float64,1}
    Tper::Float64
    Mup::Array{Float64,1}
    theta_av::Float64
    theta_sd::Float64
    whisk_amp_avr::Float64
    whisk_setup::Float64
    
    avr_val() = new()
end

mutable struct net_par
    Npop::Int64
    initial_val::Array{Float64,1}

    I0_min::Float64
    I0_max::Float64
    I0_threshold::Float64
    I0::Float64
    taua::Float64
    taus::Float64
    beta::Float64
    gamma::Float64
    gKm::Float64
    DelV::Float64
    ginter::Float64
    gintra::Float64

    DelV_F::Float64
    gFr::Float64

    inhibition_factor::Float64 

    I0_threshold_F::Float64
    I0_F::Float64
    taua_F::Float64
    beta_F::Float64
    gamma_F::Float64
    gKm_F::Float64
 
    grB::Float64
    IrB::Float64
    TB::Float64
    Tup::Float64

    Ja::Float64
    Jintra::Float64
    Jinter::Float64

    Ja_F::Float64
    JFr::Float64

    tauw::Float64
    
    net_par() = new()
end

mutable struct run_par
    epsilon::Float64
    Tall::Float64
    deltat::Float64
    NT::Int64
    method::Char
    tstat::Float64
    twrite::Int64
    tmcol::Float64
    incond::Char
    smforce::Char
    sm::Bool
    n_hist_store::Int64
    open_file_for_par::Bool
    min_theta_for_osci::Float64
    
#   for tvn.jl
    ncut::Int64
    eps_from_bif::Float64
    parmax::Float64
    
    run_par() = new()
end

mutable struct pb_str
    Tphase::Float64
    s_preBot::Float64
    s_preBot_o::Float64
    up::Bool
    n_preBot_cycle::Int64
    m_dat_cycle::Int64
    dat_cycle::Array{Array{Float64,1},1}

    pb_str() = new()
end

mutable struct stat_str
    hist_store_ar::Array{Float64,1}
    peak_time_ar::Array{Float64,1}
    Mav::Array{Float64,1}
    Mavt::Array{Float64,1}
    Msd::Array{Float64,1}

    theta_cb::CircularBuffer{Float64}
    theta_min_ar::Array{Array{Float64,1},1}
    theta_max_ar::Array{Array{Float64,1},1}

    theta_sum::Float64
    theta_sumt::Float64

    stat_str() = new()
end

mutable struct fl_st
    col::IOStream
    grb::IOStream
    pek::IOStream
    avr::IOStream
    acl::IOStream
    the::IOStream
    out::IOStream
    wsk::IOStream
    tim::IOStream
    
    fl_st() = new()
end
