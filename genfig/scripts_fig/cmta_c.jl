# This program plots the bursting time period Tper, Mfr, CV, CV2 and whisking
# amplitude and baseline computed by simulations and analytically when possible.

using Formatting
using Printf

mutable struct suffix_str
    s::String
    c::String
    t::String
    
    suffix_str() = new()
end

mutable struct file_names
    tper_s::String
    wsk::String
    tper_ct::String
    mfr_s::String
    mfr_ct::String
    cvt::String
    mfr_t::String
    mfr_t_scale::String
    wsk_a::String
    wsk_amp_ct::String
    wsk_set_ct::String
    wsk_amp_t::String
    wsk_amp_t_scale::String
    wsk_set_t::String
    wsk_set_t_scale::String
    arrows::String
    arrow_cvt::String
    
    file_names() = new()
end

mutable struct fl_st
    dats::IOStream
    datsw::IOStream
    datc::IOStream
    datt::IOStream
    tper_s::IOStream
    wsk::IOStream
    tper_ct::IOStream
    mfr_s::IOStream
    mfr_ct::IOStream
    cvt::IOStream
    mfr_t::IOStream
    mfr_t_scale::IOStream
    wsk_a::IOStream
    wsk_amp_ct::IOStream
    wsk_set_ct::IOStream
    wsk_amp_t::IOStream
    wsk_amp_t_scale::IOStream
    wsk_set_t::IOStream
    wsk_set_t_scale::IOStream
    arrows::IOStream
    arrow_cvt::IOStream
    
    fl_st() = new()
end

mutable struct sim_ar
    Tper::Array{Array{Float64,1},1}    
    Tper_sim_mean_sd::Array{Array{Float64,1},1}
    Mfr::Array{Array{Array{Float64,1},1},1}
    Mfr_sim_mean_sd::Array{Array{Array{Float64,1},1},1}
    CV::Array{Array{Array{Float64,1},1},1}
    CV_sim_mean_sd::Array{Array{Array{Float64,1},1},1}
    CVt::Array{Array{Array{Float64,1},1},1}
    CVt_sim_mean_sd::Array{Array{Array{Float64,1},1},1}
    wsk_set::Array{Array{Array{Float64,1},1},1}
    wsk_sim_set_mean_sd::Array{Array{Array{Float64,1},1},1}
    wsk_amp::Array{Array{Array{Float64,1},1},1}
    wsk_sim_amp_mean_sd::Array{Array{Array{Float64,1},1},1}
    wsk_a_set::Array{Array{Array{Float64,1},1},1}
    wsk_a_sim_set_mean_sd::Array{Array{Array{Float64,1},1},1}
    wsk_a_amp::Array{Array{Array{Float64,1},1},1}
    wsk_a_sim_amp_mean_sd::Array{Array{Array{Float64,1},1},1}
    
    sim_ar() = new()
end

# This program opens the input and output files.
function open_files_one_parameter(suffixstr, filenames, fl)

    file_dat_s = suffixstr.s * "/irt.avr"
    println("file_dat_s=", file_dat_s)
    fl.dats = open(file_dat_s, "r")

    file_dat_sw = suffixstr.s * "/irt.aws"
    println("file_dat_sw=", file_dat_sw)
    fl.datsw = open(file_dat_sw, "r")

    file_dat_t = "../../rate_num/datfig/$(suffixstr.c)/tvnc.res"
#   file_dat_t = "../../rate_ode/tv.res." * suffixstr.t
    println("file_dat_t=", file_dat_t)
    fl.datt = open(file_dat_t, "r")

    file_dat_c = "../../rate_num/datfig/$(suffixstr.c)/rf.avr"
#   file_dat_c = "../../rate_num/datfig/$(suffixstr.c)/numr.avr"
#   file_dat_c = "../../rate_num/dat_rf/b1/rf.avr"
    println("file_dat_c=", file_dat_c)
    fl.datc = open(file_dat_c, "r")

    filenames.wsk_amp_t = "../../rate_num/datfig/$(suffixstr.c)/tvnc.amp"
    filenames.wsk_amp_t_scale =
      "../../rate_num/datfig/$(suffixstr.c)/tvnc.amp.scale"
    println("file_dat_amp_t=", filenames.wsk_amp_t)
    fl.wsk_amp_t = open(filenames.wsk_amp_t, "r")
    fl.wsk_amp_t_scale = open(filenames.wsk_amp_t_scale, "w")

    filenames.wsk_set_t = "../../rate_num/datfig/$(suffixstr.c)/tvnc.set"
    filenames.wsk_set_t_scale =
      "../../rate_num/datfig/$(suffixstr.c)/tvnc.set.scale"
    println("file_dat_set_t=", filenames.wsk_set_t)
    fl.wsk_set_t = open(filenames.wsk_set_t, "r")
    fl.wsk_set_t_scale = open(filenames.wsk_set_t_scale, "w")

    filenames.tper_s  = suffixstr.s * "/cmt.tper_s"
    filenames.wsk = suffixstr.s * "/cmt.wsk"
    filenames.tper_ct = suffixstr.s * "/cmt.tper_ct"
    filenames.mfr_s  = suffixstr.s * "/cmt.mfr_s"
    filenames.mfr_ct = suffixstr.s * "/cmt.mfr_ct"
#   filenames.mfr_t = "../../rate_ode/tv.mfr." * suffixstr.t
    filenames.mfr_t = "../../rate_num/datfig/$(suffixstr.c)/tvnc.mfr"
    filenames.mfr_t_scale =
      "../../rate_num/datfig/$(suffixstr.c)/tvnc.mfr.scale"
    filenames.cvt    = suffixstr.s * "/cmt.cvt"
    filenames.wsk_a = suffixstr.s * "/cmt.wsk_a"
    filenames.wsk_amp_ct = suffixstr.s * "/cmt.wsk_amp_ct"
    filenames.wsk_set_ct = suffixstr.s * "/cmt.wsk_set_ct"
    filenames.arrows = suffixstr.s * "/arrows.xx"
    filenames.arrow_cvt = suffixstr.s * "/arrow_cvt.xx"

    fl.tper_s  = open(filenames.tper_s , "w")
    fl.wsk = open(filenames.wsk, "w")
    fl.tper_ct = open(filenames.tper_ct, "w")
    fl.mfr_s  = open(filenames.mfr_s , "w")
    fl.mfr_ct = open(filenames.mfr_ct, "w")
    fl.mfr_t = open(filenames.mfr_t, "r")
    fl.mfr_t_scale = open(filenames.mfr_t_scale, "w")    
    fl.cvt    = open(filenames.cvt,    "w")
    fl.wsk_set_t = open(filenames.wsk_set_t, "r")
    fl.wsk_set_t_scale = open(filenames.wsk_set_t_scale, "w")
    fl.wsk_a = open(filenames.wsk_a, "w")
    fl.wsk_amp_ct = open(filenames.wsk_amp_ct, "w")
    fl.wsk_set_ct = open(filenames.wsk_set_ct, "w")
    fl.arrows = open(filenames.arrows, "w")
    fl.arrow_cvt = open(filenames.arrow_cvt, "w")

end #open_files_one_parameter

# This program closes the input and output files.
function close_files_one_parameter(fl)

    close(fl.dats)
    close(fl.datsw)
    close(fl.datc)
    close(fl.datt)
    close(fl.tper_s)
    close(fl.wsk)
    close(fl.tper_ct)
    close(fl.mfr_s)
    close(fl.mfr_ct)
    close(fl.cvt)
    close(fl.wsk_a)
    close(fl.wsk_amp_ct)
    close(fl.wsk_set_ct)
    close(fl.arrows)
    close(fl.arrow_cvt)
    
end #close_files_one_parameter

# This function reads the data file fread.
function read_data_file(col_ar, fread)

    dat_ar = []
    for line in eachline(fread)
        if occursin(r"[0-9]", line)
            line_split = split(line)
	    one_record = []
	    for index in col_ar
	        if index <= length(line_split)
	            if occursin(".", line_split[index])
                        push!(one_record, parse(Float64, line_split[index]))
                    else
	                push!(one_record, parse(Int64, line_split[index]))
	            end
		else
		    push!(one_record, 0.0)
		end
	    end
	    push!(dat_ar, one_record)
	end
    end

    return dat_ar
end #read_data_file

# This function reads the data files into data arrays
function read_several_data_files(fl)

#   1_1:par
#   2_11:Mfr_1 3_12:CV_1 4_13:CV2_1
#   5_21:Mfr_2 6_22:CV_2 7_23:CV2_2
#   8_32:mode 9_34:Tper
    dats_ar = read_data_file((1, 11, 12, 13, 21, 22, 23, 32, 34, 9, 19),
      fl.dats)
#   seekstart(fl.dats)

    datsw_ar = read_data_file((1, 4, 8, 10, 11), fl.datsw) # 8 instead of 5

    datt_ar = read_data_file((1, 2), fl.datt)
    [record[2] *= 2.0 for record in datt_ar]
    datc_ar = read_data_file((1, 4, 6, 10, 15, 16), fl.datc)

    return dats_ar, datsw_ar, datt_ar, datc_ar
end #read_several_data_files

# This function returns an array pointing to the beginning and end of records
# with the same parameter.
function find_eq_par_ar(eq_par_ar, dats_ar)

#   find array of similar parameters.

    irecord = 0
    eq_par = [1]
    while irecord < length(dats_ar)
        irecord +=1
	if abs(dats_ar[irecord][1] - dats_ar[eq_par[1]][1]) > epsilon
	    push!(eq_par, irecord-1)
	    push!(eq_par_ar, eq_par)
	    eq_par = [irecord]
	end
    end
    push!(eq_par, length(dats_ar))
    push!(eq_par_ar, eq_par)
#   println("eq_par_ar=", eq_par_ar)

end #find_eq_par_ar

# This function defines simar.
function define_sim_ar(simar)

   simar.Tper = Array{Array{Float64,1},1}(undef,0)
   simar.Tper_sim_mean_sd = Array{Array{Float64,1},1}(undef,0)
   
   simar.Mfr = Array{Array{Array{Float64,1},1},1}(undef,2)
   [simar.Mfr[ipop] = Array{Array{Float64,1},1}(undef,0) for ipop in 1:2]
   simar.Mfr_sim_mean_sd = Array{Array{Array{Float64,1},1},1}(undef,2)
   [simar.Mfr_sim_mean_sd[ipop] = Array{Array{Float64,1},1}(undef,0) for
     ipop in 1:2]
     
   simar.CV = Array{Array{Array{Float64,1},1},1}(undef,2)
   [simar.CV[ipop] = Array{Array{Float64,1},1}(undef,0) for ipop in 1:2]
   simar.CV_sim_mean_sd = Array{Array{Array{Float64,1},1},1}(undef,2)
   [simar.CV_sim_mean_sd[ipop] = Array{Array{Float64,1},1}(undef,0) for
     ipop in 1:2]

   simar.CVt = Array{Array{Array{Float64,1},1},1}(undef,2)
   [simar.CVt[ipop] = Array{Array{Float64,1},1}(undef,0) for ipop in 1:2]
   simar.CVt_sim_mean_sd = Array{Array{Array{Float64,1},1},1}(undef,2)
   [simar.CVt_sim_mean_sd[ipop] = Array{Array{Float64,1},1}(undef,0) for
     ipop in 1:2]
   
   simar.wsk_set = Array{Array{Array{Float64,1},1},1}(undef,2)
   [simar.wsk_set[ipop] = Array{Array{Float64,1},1}(undef,0) for ipop in 1:2]
   simar.wsk_sim_set_mean_sd = Array{Array{Array{Float64,1},1},1}(undef,2)
   [simar.wsk_sim_set_mean_sd[ipop] = Array{Array{Float64,1},1}(undef,0) for
     ipop in 1:2]
     
   simar.wsk_amp = Array{Array{Array{Float64,1},1},1}(undef,2)
   [simar.wsk_amp[ipop] = Array{Array{Float64,1},1}(undef,0) for ipop in 1:2]
   simar.wsk_sim_amp_mean_sd = Array{Array{Array{Float64,1},1},1}(undef,2)
   [simar.wsk_sim_amp_mean_sd[ipop] = Array{Array{Float64,1},1}(undef,0) for
     ipop in 1:2]
     
   simar.wsk_a_set = Array{Array{Array{Float64,1},1},1}(undef,2)
   [simar.wsk_a_set[ipop] = Array{Array{Float64,1},1}(undef,0) for ipop in 1:2]
   simar.wsk_a_sim_set_mean_sd = Array{Array{Array{Float64,1},1},1}(undef,2)
   [simar.wsk_a_sim_set_mean_sd[ipop] = Array{Array{Float64,1},1}(undef,0) for
     ipop in 1:2]
     
   simar.wsk_a_amp = Array{Array{Array{Float64,1},1},1}(undef,2)
   [simar.wsk_a_amp[ipop] = Array{Array{Float64,1},1}(undef,0) for ipop in 1:2]
   simar.wsk_a_sim_amp_mean_sd = Array{Array{Array{Float64,1},1},1}(undef,2)
   [simar.wsk_a_sim_amp_mean_sd[ipop] = Array{Array{Float64,1},1}(undef,0) for
     ipop in 1:2]
     
end #define_sim_ar

# This function computes the average and standard deviation for repetitions for
# one parameter.
function compute_mean_sd_one_par(data_array, icol, jline_min, jline_max)

#   println("jline_min=", jline_min, " jline_max=", jline_max)

    sum = 0.0
    sum2 = 0.0
    for jl in jline_min:jline_max
        xx = data_array[jl][icol]
        sum += xx
        sum2 += xx * xx
    end

    mean = sum  / (1 + jline_max - jline_min)
    diff = sum2 / (1 + jline_max - jline_min) - mean * mean
    
    if jline_max == jline_min
        sd = 0.0
    elseif diff >= 0.0
        sd = sqrt(diff)
    elseif diff > -epsilon
        sd = 0.0
    else
        sd = -999.9
    end
    
    return mean, sd
end #compute_mean_sd_one_par

# This function computes the average and standard deviation of one data column.
function compute_mean_sd_one_col(mean_sd_ar, data_array, icol)

    icol_par = 1
    len_data = length(data_array)
    jline_min = 1
    jline_max = 1
    while jline_max < len_data
        diff = abs(data_array[jline_max+1][icol_par] -
	           data_array[jline_max  ][icol_par])
        if diff < epsilon
	    jline_max += 1
        else
	    mean, sd = compute_mean_sd_one_par(data_array, icol,jline_min,
	      jline_max)
	    push!(mean_sd_ar, [data_array[jline_max  ][icol_par], mean, sd])
	    jline_min = jline_max + 1
	    jline_max = jline_min
	end
    end

    if length(data_array) > 0
        mean, sd = compute_mean_sd_one_par(data_array, icol,jline_min,
	  jline_max)
         push!(mean_sd_ar, [data_array[jline_max  ][icol_par], mean, sd])
    end

end #compute_mean_sd_one_col

# This function extracts Tper, Mfr, CV and CV2 for the spiking neuron model.
# For Tper, it considers only parameters with majority of records that yield
# bursting, mode=2.
function find_T_M_CV_sim(simar, eq_par_simar, dats_ar, datsw_ar, x_scale, fl)

#   Tper
#   Parameters are considered to belong to the periodic regime if the number
#   of periodic redord is larger than ratio_per_min-epsilon
    ratio_per_min = 0.6
    nmode = 8
    nTper = 9
    nMfr1 = 2
    nMfr2 = 5
    nCV1 = 3
    nCV2 = 6
    nCVt1 = 4
    nCVt2 = 7
    nCVt1_no_bur = 10
    nCVt2_no_bur = 11

    for eq_par in eq_par_simar
        periodic_ar = []
        for irecord in eq_par[1]:eq_par[2]
	    record =  dats_ar[irecord]
#           println("record=", record)
	    if record[nmode] == 2 &&
	      record[nTper] > 0.5 * 4.0 * 1.0 / (record[nMfr1] + record[nMfr2])
	        push!(periodic_ar, irecord)
	    end
	end
	ratio_per = 1.0 * length(periodic_ar) / (eq_par[2] - eq_par[1] + 1) 

        if ratio_per > ratio_per_min - epsilon
	    for irecord in periodic_ar
	        record = dats_ar[irecord]
#		println("irecord=", irecord, " record=", record)
	        push!(simar.Tper, [record[1], record[nTper]])
	    end
	end
    end

#   Averaging Tper over realizations with the same parameters.
    if length(simar.Tper) > 0
        compute_mean_sd_one_col(simar.Tper_sim_mean_sd, simar.Tper, 2)
	for record in simar.Tper_sim_mean_sd
	    if record[2] > 0.0
	        ff = 1000.0 / record[2]
		Delff = 1000.0 * record[3] / (record[2]^2)
	        printfmtln(fl.tper_s, "{1:f} {2:f} {3:f}", x_scale * record[1], ff, Delff)
#		  1000.0 /record[2], record[3])
	    end
	end
    end

#   Averaging Mfr over realizations with the same parameters.
    fr_considered_zero = 0.005
    for (record, recordw) in zip(dats_ar, datsw_ar)
        if record[nMfr1] <= fr_considered_zero &&
	  record[nMfr2] <= fr_considered_zero
	    CVt1 = -997.7
	    CVt2 = -997.7
	elseif record[nMfr1] <= fr_considered_zero
	    CVt1 = -997.7
	    CVt2 = record[nCVt2_no_bur]
	elseif record[nMfr2] <= fr_considered_zero
	    CVt1 = record[nCVt1_no_bur]
	    CVt2 = -997.7
	else
##	    printfmt("r1={1:f} T={2:f} M={3:f} {4:f} 1/M={5:f} {6:f}",
##	      x_scale * record[1], record[nTper], record[nMfr1], record[nMfr2],
##	      1.0/record[nMfr1], 1.0/record[nMfr2])
##	    printfmt(" CVt1={1:f} {2:f} CVt2={3:f} {4:f}", record[nCVt1],
##	      record[nCVt1_no_bur], record[nCVt2], record[nCVt2_no_bur])
	    if record[nTper] > 2.0 * 1.0 / record[nMfr1]
	        CVt1 = record[nCVt1]
            else
	        CVt1 = record[nCVt1_no_bur]
            end
	    if record[nTper] > 2.0 * 1.0 / record[nMfr2]
	        CVt2 = record[nCVt2]
            else
	        CVt2 = record[nCVt2_no_bur]
            end
	end
	
##	printfmtln(" A CVt1={1:f} CVt2={2:f}", CVt1, CVt2)	

	if record[nMfr1] >= record[nMfr2]
	    push!(simar.Mfr[1], [record[1], record[nMfr1]])
	    push!(simar.Mfr[2], [record[1], record[nMfr2]])
	    push!(simar.CV[1], [record[1], record[nCV1]])
	    push!(simar.CV[2], [record[1], record[nCV2]])
#           record[nCVt1]])
	    CVt1 > 0.0 ? push!(simar.CVt[1], [record[1], CVt1]) : nothing
#           # record[nCVt2]])
	    CVt2 > 0.0 ? push!(simar.CVt[2], [record[1], CVt2]) : nothing
	    push!(simar.wsk_set[1], [recordw[1], recordw[2]])
	    push!(simar.wsk_amp[1], [recordw[1], recordw[3]])
	    push!(simar.wsk_a_set[1], [recordw[1], recordw[4]])
	    push!(simar.wsk_a_amp[1], [recordw[1], recordw[5]])
	else
	    push!(simar.Mfr[1], [record[1], record[nMfr2]])
	    push!(simar.Mfr[2], [record[1], record[nMfr1]])
	    push!(simar.CV[1], [record[1], record[nCV2]])
	    push!(simar.CV[2], [record[1], record[nCV1]])
#           record[nCVt2]])
            CVt2 > 0.0 ? push!(simar.CVt[1], [record[1], CVt2]) : nothing
#           record[nCVt1]])
	    CVt1 > 0.0 ? push!(simar.CVt[2], [record[1], CVt1]) : nothing
	    push!(simar.wsk_set[2], [recordw[1], recordw[2]])
	    push!(simar.wsk_amp[2], [recordw[1], recordw[3]])
	    push!(simar.wsk_a_set[2], [recordw[1], recordw[4]])
	    push!(simar.wsk_a_amp[2], [recordw[1], recordw[5]])
         end
	 
#	 if abs(record[1] - 15.0) < 0.001
#	     println(record)
#	 end
    end

#   Averaging Mfr over realizations with the same parameters.
    for ipop in 1:2
        compute_mean_sd_one_col(simar.Mfr_sim_mean_sd[ipop], simar.Mfr[ipop], 2)
	for record in simar.Mfr_sim_mean_sd[ipop]
	    printfmtln(fl.mfr_s, "{1:f} {2:f} {3:f}", x_scale * record[1],
	      1000.0*record[2], 1000.0*record[3])
	end
	println(fl.mfr_s, "  ")
    end
    
#   Averaging CV over realizations with the same parameters.
    CV_print = false
    if CV_print
        for ipop in 1:2
            compute_mean_sd_one_col(simar.CV_sim_mean_sd[ipop], simar.CV[ipop],
	      2)
	    for (record_Mfr, record_CV) in
	      zip(simar.Mfr_sim_mean_sd[ipop], simar.CV_sim_mean_sd[ipop])
	        if record_Mfr[2] > 0.005
	            printfmtln(fl.cvt, "{1:f} {2:f} {3:f}", x_scale * record_CV[1],
		      record_CV[2], record_CV[3])
                else
	            printfmtln(fl.cvt, "{1:f} -1.0 0.0", x_scale * record_CV[1])
                end
	    end
	    println(fl.cvt, "  ")
        end
    end

#   Averaging CVt over realizations with the same parameters.
    for ipop in 1:2
        compute_mean_sd_one_col(simar.CVt_sim_mean_sd[ipop], simar.CVt[ipop], 2)
	for (record_Mfr, record_CVt) in
	  zip(simar.Mfr_sim_mean_sd[ipop], simar.CVt_sim_mean_sd[ipop])
	    if record_Mfr[2] > 0.005
	        printfmtln(fl.cvt, "{1:f} {2:f} {3:f}", x_scale * record_CVt[1],
		  record_CVt[2], record_CVt[3])
            else
	        printfmtln(fl.cvt, "{1:f} -1.0 0.0", x_scale * record_CVt[1])
            end
	end
	println(fl.cvt, "  ")
     end

    angle_factor = 0.12
#   Averaging wsk_set over realizations with the same parameters.
    for ipop in 1:2
        compute_mean_sd_one_col(simar.wsk_sim_set_mean_sd[ipop],
	  simar.wsk_set[ipop], 2)
#	println("ipop=", ipop)
#	println("wsk_set=", simar.wsk_set[ipop])
#	println("mean=", simar.wsk_sim_set_mean_sd[ipop])
	for (record_Mfr, record_wsk_set) in
	  zip(simar.Mfr_sim_mean_sd[ipop], simar.wsk_sim_set_mean_sd[ipop])
	    printfmtln(fl.wsk, "{1:f} {2:f} {3:f}", x_scale * record_wsk_set[1],
              angle_factor * record_wsk_set[2],
	      angle_factor * record_wsk_set[3])
	end
	    println(fl.wsk, "  ")
     end

#   Averaging wsk_amp over realizations with the same parameters.
    for ipop in 1:2
        compute_mean_sd_one_col(simar.wsk_sim_amp_mean_sd[ipop],
	  simar.wsk_amp[ipop], 2)
	for (record_Mfr, record_wsk_amp) in
	  zip(simar.Mfr_sim_mean_sd[ipop], simar.wsk_sim_amp_mean_sd[ipop])
	    printfmtln(fl.wsk, "{1:f} {2:f} {3:f}", x_scale * record_wsk_amp[1],
              angle_factor * record_wsk_amp[2],
	      angle_factor * record_wsk_amp[3])
	end
	    println(fl.wsk, "  ")
     end

#   Averaging wsk_a_set over realizations with the same parameters.
    for ipop in 1:2
        compute_mean_sd_one_col(simar.wsk_a_sim_set_mean_sd[ipop],
	  simar.wsk_a_set[ipop], 2)
#	println("ipop=", ipop)
#	println("wsk_a_set=", simar.wsk_a_set[ipop])
#	println("mean=", simar.wsk_a_sim_set_mean_sd[ipop])
	for (record_Mfr, record_wsk_a_set) in
	  zip(simar.Mfr_sim_mean_sd[ipop], simar.wsk_a_sim_set_mean_sd[ipop])
	    printfmtln(fl.wsk_a, "{1:f} {2:f} {3:f}", x_scale * record_wsk_a_set[1],
              angle_factor * record_wsk_a_set[2],
	      angle_factor * record_wsk_a_set[3])
	end
	    println(fl.wsk_a, "  ")
     end

#   Averaging wsk_a_amp over realizations with the same parameters.
    for ipop in 1:2
        compute_mean_sd_one_col(simar.wsk_a_sim_amp_mean_sd[ipop],
	  simar.wsk_a_amp[ipop], 2)
	for (record_Mfr, record_wsk_a_amp) in
	  zip(simar.Mfr_sim_mean_sd[ipop], simar.wsk_a_sim_amp_mean_sd[ipop])
	    printfmtln(fl.wsk_a, "{1:f} {2:f} {3:f}", x_scale * record_wsk_a_amp[1],
              0.5 * angle_factor * record_wsk_a_amp[2],
	      0.5 * angle_factor * record_wsk_a_amp[3])
	end
	    println(fl.wsk_a, "  ")
     end

end #find_T_M_CV_sim

# This function computes the average of whisking setpoint and amplitude
# over realizations with the same parameter for the spiking neuron model.
#function find_W_sim(simar, datsw_ar, fl)
#
#   for iwsk_var in 1:2
#       compute_mean_sd_one_col(simar.wsk_sim_set_mean_sd[iwsk_var], datsw_ar,
#	  iwsk_var+1)
#	for record in simar.wsk_sim_set_mean_sd[iwsk_var]
#	    printfmtln(fl.wsk, "{1:f} {2:f} {3:f}", x_scale * record[1], record[2],
#	      record[4]) # record[3]: <sigma_theta>, record[4]: amp
#	end
#	println(fl.wsk, "  ")
#   end

#end #find_W_sim(datsw_ar, fl)

# This function writes the data computed numerically by simulating the rate
# model.
function write_rate_num(datc_ar, x_scale, fl)

    iTper = 4
    for record in datc_ar
        if record[iTper] > 50.0 # 0.0
	    printfmtln(fl.tper_ct, "{1:f} {2:f}", x_scale * record[1],
	      1000.0 / record[iTper])
	end
    end
    
    for (index, record) in enumerate(datc_ar)
        if abs(record[2] - record[3]) <= 0.01
            println(fl.mfr_ct, x_scale * record[1], " ",
	      1000.0 * 0.5 * (record[2] + record[3]))
	else
            println(fl.mfr_ct, x_scale * record[1], " ",
	      1000.0 * max(record[2], record[3]))
	end
    end
    
    println(fl.mfr_ct, "  ")
    
    for (index, record) in enumerate(datc_ar)
        if abs(record[2] - record[3]) > 0.01
            println(fl.mfr_ct, x_scale * record[1], " ",
	      1000.0 * min(record[2], record[3]))
	else
	    if index < length(datc_ar)
	        record_next = datc_ar[index+1]
		if abs(record_next[2] - record_next[3]) > 0.01
                    println(fl.mfr_ct, x_scale * record[1], " ",
	              1000.0 * 0.5 * (record[2] + record[3]))
		end
	    end
	end
    end

    i_wsk_amp = 5
    for record in datc_ar
	printfmtln(fl.wsk_amp_ct, "{1:f} {2:f}", x_scale * record[1],
	  0.5 * record[i_wsk_amp])
    end
    
    i_wsk_set = 6
    for (irecord, record) in enumerate(datc_ar)
	printfmtln(fl.wsk_set_ct, "{1:f} {2:f}", x_scale * record[1], record[i_wsk_set])
	if irecord < length(datc_ar)
	    if record[i_wsk_set] > datc_ar[irecord+1][i_wsk_set]
	        println(fl.wsk_set_ct, "  ")	        
	    end
	end
    end

    for irecord in 1:length(datc_ar)-1
        if datc_ar[irecord][i_wsk_amp] > 10.0 &&
	  datc_ar[irecord+1][i_wsk_amp] < 1.0e-10
	      println(fl.wsk_set_ct, "  ")
#	      printfmtln(fl.wsk_set_ct, "{1:f} {2:f}", datc_ar[irecord+1][1],
#	        datc_ar[irecord+1][i_wsk_set])
	      printfmtln(fl.wsk_set_ct, "{1:f} 0.0", datc_ar[irecord+1][1])
              printfmtln(fl.wsk_set_ct, "{1:f} 0.0", datc_ar[end][1])
	end
    end
    
end #write_rate

function write_rate_th(datt_ar, x_scale, fl)

    println(fl.tper_ct, "  ")
    for record in datt_ar
        printfmtln(fl.tper_ct, "{1:f} {2:f}", x_scale * record[1], 1000.0 / record[2])
    end

end #write_rate_th

# This function writes the arrows data in the arrows file.
function write_arrows(file_num, fl)

     if file_num == 1
         println(fl.arrows, "6.0 0.0")
	 println(fl.arrows, "6.0 0.6")
	 println(fl.arrows, "  ")
	 println(fl.arrows, "6.0 0.6")
	 println(fl.arrows, "  ")
         println(fl.arrows, "12.0 0.0")
	 println(fl.arrows, "12.0 0.6")
	 println(fl.arrows, "  ")
	 println(fl.arrows, "12.0 0.6")
	 println(fl.arrows, "  ")
         printfmtln(fl.arrows, "{1:f} 0.0", 3.148666 / 25.0)
	 printfmtln(fl.arrows, "{1:f} 0.6", 3.148666 / 25.0)
	 printfmtln(fl.arrows, "{1:f} 0.6", 8.527547 / 25.0)
         printfmtln(fl.arrows, "{1:f} 0.0", 8.527547 / 25.0)
         printfmtln(fl.arrows, "{1:f} 0.0", 3.148666 / 25.0)
	 println(fl.arrows, "  ")
	 
     elseif file_num == 2
         println(fl.arrows, "  ")
     end
     
     println(fl.arrow_cvt, "6.0 0.0")
     println(fl.arrow_cvt, "5.6 0.05")
     println(fl.arrow_cvt, "6.0 0.04")
     println(fl.arrow_cvt, "6.4 0.05")
     println(fl.arrow_cvt, "6.0 0.0")

end #write_arrows

# This function scales the theoretical results by multiplying the x-axis by
# x_scale.
function scale_file(fread, fwrite, x_scale)

    for line in eachline(fread)
        if occursin(r"([0-9])", line)
	    line_split_float = parse.(Float64, split(line))
	    printfmtln(fwrite, "{1:f} {2:f}", x_scale * line_split_float[1],
	      line_split_float[2])
	else
	    println(fwrite, line)
	end
    end

    close(fwrite)
end

# This function process the data from a set of simulations in which
# one parameter is podified.
function process_modify_one_parameter(suffixstr, filenames, fl, file_num,
  x_scale)

    open_files_one_parameter(suffixstr, filenames, fl)

#   Reading the data files into data arrays
    dats_ar, datsw_ar, datt_ar, datc_ar =
      read_several_data_files(fl)

#   Spiking neuron model, vIRT
    eq_par_simar = Array{Array{Int64,1}}(undef,0)
    find_eq_par_ar(eq_par_simar, dats_ar)
    simar = sim_ar()
    define_sim_ar(simar)
    find_T_M_CV_sim(simar, eq_par_simar, dats_ar, datsw_ar, x_scale, fl)

#   Spiking neuron model, whisker
#   find_W_sim(simar, datsw_ar, fl)

#   Rate model, simulations
    write_rate_num(datc_ar, x_scale, fl)

#   Rate model, theory
    write_rate_th(datt_ar, x_scale, fl)
    scale_file(fl.mfr_t, fl.mfr_t_scale, x_scale)
    scale_file(fl.wsk_amp_t, fl.wsk_amp_t_scale, x_scale)
    scale_file(fl.wsk_set_t, fl.wsk_set_t_scale, x_scale)
    write_arrows(file_num, fl)

    close_files_one_parameter(fl)

end #process_modify_one_parameter

# This function plots a column of graphs for one parameter set.
function plot_one_parameter_set(xm_cmd, filenames, ngraph)

    mgraph = 5
    igraph = mgraph * ngraph
    
    push!(xm_cmd, "-graph", string(igraph))
    push!(xm_cmd, filenames.arrows )

    igraph += 1
    if ngraph == 0
        push!(xm_cmd, "-graph", string(igraph))
        push!(xm_cmd, filenames.arrows )
    end
    
    igraph += 1
    push!(xm_cmd, "-graph", string(igraph))
    push!(xm_cmd, "-settype", "xydx", filenames.mfr_s )
    push!(xm_cmd, "-settype", "xy",   filenames.mfr_ct)
    push!(xm_cmd, "-settype", "xy",   filenames.mfr_t_scale )

    igraph += 1
    push!(xm_cmd, "-graph", string(igraph))
    push!(xm_cmd,  "-settype", "xydx", filenames.tper_s )
    push!(xm_cmd,  "-settype", "xy",   filenames.tper_ct)
#   "-settype", "xydy",

    igraph += 1
    push!(xm_cmd, "-graph", string(igraph))
    push!(xm_cmd, "-settype", "xydx", filenames.wsk_a)
    push!(xm_cmd, "-settype", "xy", filenames.wsk_amp_ct)
    push!(xm_cmd, "-settype", "xy", filenames.wsk_amp_t_scale)
    
    igraph += 1
    push!(xm_cmd, "-graph", string(igraph))
    push!(xm_cmd, "-settype", "xydx", filenames.wsk_a)
    push!(xm_cmd, "-settype", "xy", filenames.wsk_set_ct)
    push!(xm_cmd, "-settype", "xy", filenames.wsk_set_t_scale)
    
end #plot_one_parameter_set

# This function plots a column of graphs for one parameter set.
function plot_one_parameter_set_a(xm_cmd, filenames, ngraph)

    mgraph = 1
    igraph = mgraph * ngraph
    
    if ngraph == 0
        push!(xm_cmd, "-graph", string(igraph))
        push!(xm_cmd, filenames.arrows )
    end
    
    igraph += 1
    push!(xm_cmd, "-graph", string(igraph))
    push!(xm_cmd, "-settype", "xydy", filenames.cvt)
    if ngraph == 0
        push!(xm_cmd, "-settype", "xy", filenames.arrow_cvt)
    end
    
#   igraph += 1
#   push!(xm_cmd, "-graph", string(igraph))
#   push!(xm_cmd, "-settype", "xydx", filenames.wsk)
    
#   igraph += 1
#   push!(xm_cmd, "-graph", string(igraph))
#   push!(xm_cmd, "-settype", "xydx", filenames.wsk)
    
end #plot_one_parameter_set_a

# This function plots the file fcol_name
function plot_curve(suffix_s_a, filenames_a, suffix_s_b, filenames_b,
  gr_fl)

    HOME = "/home/golomb"
    xm_cmd = String[]
    push!(xm_cmd, "xmgrace")
#   push!(xm_cmd, "gracebat")

    plot_one_parameter_set(xm_cmd, filenames_a, 0)
    plot_one_parameter_set(xm_cmd, filenames_b, 1)

    push!(xm_cmd, "-hdevice", "EPS", "-p",  gr_fl)
    push!(xm_cmd, "-printfile", suffix_s_a * "/cmt_c.eps")

    println(xm_cmd)
    run(`$xm_cmd`)
    
end #plot_curve

# This function plots the file fcol_name
function plot_curve_a(suffix_s_a, filenames_a, suffix_s_b, filenames_b,
  gr_fl)

    HOME = "/home/golomb"
    xm_cmd = String[]
    push!(xm_cmd, "xmgrace")
#   push!(xm_cmd, "gracebat")

    plot_one_parameter_set_a(xm_cmd, filenames_a, 0)
    plot_one_parameter_set_a(xm_cmd, filenames_b, 1)

    push!(xm_cmd, "-hdevice", "EPS", "-p",  gr_fl)
    push!(xm_cmd, "-printfile", suffix_s_a * "/cmt_c_a.eps")

    println(xm_cmd)
    run(`$xm_cmd`)
    
end #plot_curve

#main

fl_Js_0 = fl_st()
fl_Js_12 = fl_st()

suffix_Js_0 = suffix_str()
suffix_Js_12 = suffix_str()

# Specific cases
suffix_Js_0 = suffix_str()
suffix_Js_12 = suffix_str()
filenames_Js_12 = file_names()
filenames_Js_0 = file_names()

# Js = 0
suffix_Js_0.s = "dir_qsub/no_pb_x_ginter_no_intra" # c451
suffix_Js_0.c = "no_pb_x_ginter_no_intra" #c451
suffix_Js_0.t = "no_pb_x_ginter_no_intra" #c451
# Js = 12
suffix_Js_12.s = "dir_qsub/no_pb_x_ginter_no_intra_I0" # c453
suffix_Js_12.c = "no_pb_x_ginter_no_intra_I0" #c453
suffix_Js_12.t = "no_pb_x_ginter_no_intra_I0" #c453

fout = open(suffix_Js_0.s * "/cmt.out", "w")
gr_fl = "/home/golomb/shares/docu/work_ms/brainstem/irt/genfig/scripts_fig/cmta_c.gr"
gr_fl_a = "/home/golomb/shares/docu/work_ms/brainstem/irt/genfig/scripts_fig/cmta_c_a.gr"
epsilon = 1.0e-10

process_modify_one_parameter(suffix_Js_0,  filenames_Js_0 , fl_Js_0 , 1, 0.04)
process_modify_one_parameter(suffix_Js_12, filenames_Js_12, fl_Js_12, 2, 1.0 )

plot_curve(suffix_Js_0.s, filenames_Js_0, suffix_Js_12.s, filenames_Js_12,
  gr_fl)
plot_curve_a(suffix_Js_0.s, filenames_Js_0, suffix_Js_12.s, filenames_Js_12,
  gr_fl_a)

close(fout)
