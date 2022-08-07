mutable struct peak_par
    epsilon::Float64
    deltat_spk::Float64
    deltat::Float64
    nburst::Int64
    sm::Bool

    peak_par() = new()
end

mutable struct par_to_hop
    epsilon::Float64
    scan_type::Char
    open_file_for_par::Bool
    Tall::Float64
    determine_time_prebot::Char
    Tper_syn_ext::Float64
    Tup_syn_ext::Float64
    ipop_cal_wsk::Int64
    angle_cal::Char
    whisk_cal::Bool
    nonF::Int64
    tstat::Float64
    T_diff_pre_interval::Array{Float64,1}
    frac_sd_mnmx::Float64

    
    whisker::whisker_par  
    
    par_to_hop() = new()
end

mutable struct extr_val
    t::Float64
    x::Float64
    i::Int64
    mnmx::Int64
    t_extrapolate::Float64
    
    extr_val() = new()
end

mutable struct wsk_properties
    i_prebot::Int64
    i_exter_first::Int64
    i_exter_last::Int64
    amp::Array{Float64,1}

    wsk_properties() = new()
end

mutable struct aws_val
    whisk_set::Float64
    whisk_amp::Float64
    coeff_pred1::Float64
    coeff_pred2::Float64
    amp_avr::Float64
    amp_sd::Float64
    theta_max::Float64
    theta_amp::Float64
    
    aws_val() = new()
end


#mutable struct fl_st
#    col::IOStream
#    trj::IOStream
#    wsk::IOStream
#    rst::IOStream
#    avr::IOStream
#    res::IOStream
#    hei::IOStream
#    aws::IOStream
#    out::IOStream
#    
#    fl_st() = new()
#end
