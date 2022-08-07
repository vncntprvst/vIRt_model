# This program generate files from the type irt_const.jl to generate parameter
# files for simulations.

using Formatting
using Printf

# structure for modified parameter set
mutable struct mod_par
    nrun::Int64
    parar::Array{NamedTuple{(:parmin, :parmax, :npar, :nrepeat),
      Tuple{Float64, Float64, Int64, Int64}}}
	    
    mod_par() = new()
end

# This function reads the original const file.
function read_org_file(orig_file)

    fconst_org = open(dir_orig_file * "/irt_const.jl", "r")
    orig_const_ar = readlines(fconst_org)
    close(fconst_org)
#   println(orig_const_ar)

    return orig_const_ar
end #read_org_file(dir_orig_file)

# This function composes the command used to modify and write the new parameter
# file.
function compose_write_command(sim_description, array_to_interit, run_dict_str)

    cw_command = sim_description * "_const_ar = write_target_file(" * 
      array_to_interit * ", "  * sim_description * ", \""  * sim_description *
      "\", run_dict_dict[\"" * run_dict_str * "\"])"

#   println(cw_command)

    return cw_command
end #compose_write_command

# This program modifies running parameters in sval: parmin, parmax, npar,
# nrepeat.
function modify_par(const_ar, par_value)

    replace_line("sval.parmin" , string(par_value.parmin) , const_ar)
    replace_line("sval.parmax" , string(par_value.parmax) , const_ar)
    replace_line("sval.npar"   , string(par_value.npar)   , const_ar)
    replace_line("sval.nrepeat", string(par_value.nrepeat), const_ar)

    println("par_value=", par_value, " l7=", const_ar[7])
    
end #modify_par

# This function writes the target const file.
function write_target_file(original_const_ar, modify_function, dir_target_file,
  run_dict)
#  run_command, runfl)

    modpar = mod_par()
    modpar.parar = Array{NamedTuple{(:parmin, :parmax, :npar, :nrepeat),
      Tuple{Float64, Float64, Int64, Int64}}}(undef, 0)

    target_const_ar = modify_function(original_const_ar, modpar)
    
 #  generting the directory with the original sval.par...
    sub_dir_name = run_dict["dir"] * "/" * dir_target_file
    isdir(sub_dir_name) ? nothing : mkdir(sub_dir_name)
    fconst_target = open(sub_dir_name * "/irt_const.jl", "w")
    [println(fconst_target, line) for line in target_const_ar]
    close(fconst_target)
	
    println("length(parar))=", length(modpar.parar))
    if length(modpar.parar) == 0
        println(run_dict["run_fl"], run_dict["command"] * dir_target_file *
	  " &")	
    else
        for (index, par_value) in enumerate(modpar.parar)
	     dir_target_file_index = sub_dir_name * "_" *
	       string(index)
             isdir(dir_target_file_index) ? nothing :
	       mkdir(dir_target_file_index)
             fconst_target = open(dir_target_file_index * "/irt_const.jl", "w")
             modify_par(target_const_ar, par_value)
             [println(fconst_target, line) for line in target_const_ar]
             close(fconst_target)
             println(run_dict["run_fl"], run_dict["command"] *
	       dir_target_file_index)
        end
	
	if run_dict["dir"] == "dir_qsub"
	    println(num_tasks_here_fl, dir_target_file, " ", 
              length(modpar.parar))
        else
            println(num_tasks_here_qpd_fl, dir_target_file, " ", 
              length(modpar.parar))
        end
    end

    return target_const_ar
end #write_target_file

# This function replaces the value of parstr by valstr.
function replace_line(parstr, valstr, const_ar)

    iline = 1
    while iline <= length(const_ar)
        line = const_ar[iline]
        if occursin(parstr, line)
	    println("iline=", iline, " line=", line)
	    line_split = split(line, "=")
	    const_ar[iline] = line_split[1] * "= " * valstr
	    println("iline=", iline, " line=", const_ar[iline])
	    break
	end
	iline += 1
    end

    if iline > length(const_ar)
        println("iline=", iline, " > length(const_ar)=", length(const_ar))
	println("parstr=", parstr, " is not found!")
	exit(0)
    end
    
end #replace_line

# This function adds the line new_line at nline array members beyond the
# position of parstr.
function add_line(parstr, new_line, const_ar, nline)
    for iline in 1:length(const_ar)
        line = const_ar[iline]
        if occursin(parstr, line)
	    println("line=", line)
	    insert!(const_ar, iline + nline, new_line) 
	    println("line=", const_ar[iline + nline])
	    break
	end
    end
    
end #add_line

# This function deletes a line that includes parstr.
function delete_line(parstr,  const_ar)
    for iline in 1:length(const_ar)
        line = const_ar[iline]
        if occursin(parstr, line)
	    println("line to delete=", line)
            deleteat!(const_ar, iline)
	    break
	end
    end
    
end #delete_line

# This function substitutes the modified variable
function pb_v_rast_intra(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    println("  ")
    
    return target_const_ar
end #pb_v_rast_intra

# This function substitutes the modified variable
function pb_v_rast_no_intra(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line(r"^syncoupparI1I1.syn_receptor_par_ar\[1\].gsyn", "0.0",
      target_const_ar)
    replace_line(r"^syncoupparI1I2.syn_receptor_par_ar\[1\].gsyn", "0.0",
      target_const_ar)
    replace_line(r"^syncoupparI2I1.syn_receptor_par_ar\[1\].gsyn", "0.0",
      target_const_ar)
    replace_line(r"^syncoupparI2I2.syn_receptor_par_ar\[1\].gsyn", "0.0",
      target_const_ar)
    replace_line("netpar.Tper_syn_ext", "200.0", target_const_ar)
    replace_line("netpar.Trand_syn_ext", "10.0", target_const_ar)
    
    println("  ")
    
    return target_const_ar
end #pb_v_rast_no_intra

# This function substitutes the modified variable
function no_pb_v_rast_intra(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line("cellparI_mn.gsyn_ext", "0.0", target_const_ar)
    println("  ")
    
#   modpar.nrun = 0
    
    return target_const_ar
end #no_pb_v_rast_intra

# This function substitutes the modified variable
function no_pb_v_rast_no_intra(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line(r"^syncoupparI1I1.syn_receptor_par_ar\[1\].gsyn", "0.0",
      target_const_ar)
    replace_line(r"^syncoupparI1I2.syn_receptor_par_ar\[1\].gsyn", "6.0",
      target_const_ar)
    replace_line(r"^syncoupparI2I1.syn_receptor_par_ar\[1\].gsyn", "6.0",
      target_const_ar)
    replace_line(r"^syncoupparI2I2.syn_receptor_par_ar\[1\].gsyn", "0.0",
      target_const_ar)
    println("  ")
    
#   modpar.nrun = 0
    
    return target_const_ar
end #no_pb_v_rast_no_intra

# This function substitutes the modified variable
function no_pb_v_rast_no_intra(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line(r"^syncoupparI1I1.syn_receptor_par_ar\[1\].gsyn", "0.0",
      target_const_ar)
    replace_line(r"^syncoupparI1I2.syn_receptor_par_ar\[1\].gsyn", "6.0",
      target_const_ar)
    replace_line(r"^syncoupparI2I1.syn_receptor_par_ar\[1\].gsyn", "6.0",
      target_const_ar)
    replace_line(r"^syncoupparI2I2.syn_receptor_par_ar\[1\].gsyn", "0.0",
      target_const_ar)
    println("  ")
    
#   modpar.nrun = 0
    
    return target_const_ar
end #no_pb_v_rast_no_intra

# This function substitutes the modified variable
function no_pb_v_rast_intra_as(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line(r"^syncoupparI1I2.syn_receptor_par_ar\[1\].gsyn", "33.0",
      target_const_ar)
    replace_line(r"^syncoupparI2I1.syn_receptor_par_ar\[1\].gsyn", "33.0",
      target_const_ar)
    println("  ")
    
#   modpar.nrun = 0
    
    return target_const_ar
end #no_pb_v_rast_intra_as

# This function substitutes the modified variable
function no_pb_v_rast_no_intra_as(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line(r"^syncoupparI1I2.syn_receptor_par_ar\[1\].gsyn", "12.0",
      target_const_ar)
    replace_line(r"^syncoupparI2I1.syn_receptor_par_ar\[1\].gsyn", "12.0",
      target_const_ar)
    println("  ")
    
    return target_const_ar
end #no_pb_v_rast_no_intra_as

# This function substitutes the modified variable
function no_pb_x_ginter_intra(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line("sval.scan_type", "'y'", target_const_ar)
    replace_line("sval.parname",
      "[\"netpar.S_orig[1,2].syn_receptor_par_ar[1].gsyn\",",
#     "[\"syncoupparI1I2.syn_receptor_par_ar[1].gsyn\",",
      target_const_ar)
    add_line("sval.parname",
      "\"netpar.S_orig[2,1].syn_receptor_par_ar[1].gsyn = " *
      "netpar.S_orig[1,2].syn_receptor_par_ar[1].gsyn\"]",
#     "\"syncoupparI2I1.syn_receptor_par_ar[1].gsyn = " *
#     "syncoupparI1I2.syn_receptor_par_ar[1].gsyn\"]",
      target_const_ar, 1)
    replace_line("sval.parmax", "20.0", target_const_ar)
    replace_line("sval.npar", "20", target_const_ar)
    println("  ")    
    
    push!(modpar.parar, (parmin=0.0 , parmax=10.0, npar=10, nrepeat=5))
    push!(modpar.parar, (parmin=11.0, parmax=21.0, npar=10, nrepeat=5))
    push!(modpar.parar, (parmin=22.0, parmax=32.0, npar=10, nrepeat=5))
    push!(modpar.parar, (parmin=33.0, parmax=43.0, npar=10, nrepeat=5))
    
    return target_const_ar
end #no_pb_x_ginter_intra

# This function substitutes the modified variable
function no_pb_x_ginter_intra_I0(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line("sval.scan_type", "'y'", target_const_ar)
    replace_line("sval.parname",
      "[\"netpar.C_orig[1].Iapp\", \"netpar.C_orig[2].Iapp=netpar.C_orig[1].Iapp\"]",
#     "[\"cellparI1.Iapp\", \"cellparI2.Iapp=cellparI1.Iapp\"]",
      target_const_ar)
    replace_line("sval.parmax", "20.0", target_const_ar)
    replace_line("sval.npar", "20", target_const_ar)
    println("  ")
    
    push!(modpar.parar, (parmin=0.0 , parmax=10.0, npar=10, nrepeat=5))
    push!(modpar.parar, (parmin=11.0, parmax=21.0, npar=10, nrepeat=5))
    push!(modpar.parar, (parmin=22.0, parmax=32.0, npar=10, nrepeat=5))
    
    return target_const_ar
end #no_pb_x_ginter_intra_I0

# This function substitutes the modified variable
function no_pb_x_ginter_no_intra(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line("sval.scan_type", "'y'", target_const_ar)
    replace_line("sval.parname", 
      "[\"netpar.S_orig[1,2].syn_receptor_par_ar[1].gsyn\",",
#     "[\"syncoupparI1I2.syn_receptor_par_ar[1].gsyn\",",
      target_const_ar)
    add_line("sval.parname",
      "\"netpar.S_orig[2,1].syn_receptor_par_ar[1].gsyn = " *
      "netpar.S_orig[1,2].syn_receptor_par_ar[1].gsyn\"]",
#     "\"syncoupparI2I1.syn_receptor_par_ar[1].gsyn = " *
#     "syncoupparI1I2.syn_receptor_par_ar[1].gsyn\"]",
      target_const_ar, 1)
    replace_line("sval.parmax", "20.0", target_const_ar)
    replace_line("sval.npar", "20", target_const_ar)
    println("  ")
    
    push!(modpar.parar, (parmin=0.0 , parmax=10.0, npar=10, nrepeat=5))
    push!(modpar.parar, (parmin=11.0, parmax=21.0, npar=10, nrepeat=5))
     
    return target_const_ar
end #no_pb_x_ginter_no_intra

# This function substitutes the modified variable
function no_pb_x_ginter_no_intra_I0(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line("sval.scan_type", "'y'", target_const_ar)
    replace_line("sval.parname", 
      "[\"netpar.C_orig[1].Iapp\", \"netpar.C_orig[2].Iapp=netpar.C_orig[1].Iapp\"]",
#     "[\"cellparI1.Iapp\", \"cellparI2.Iapp=cellparI1.Iapp\"]",
      target_const_ar)
    replace_line("sval.parmax", "20.0", target_const_ar)
    replace_line("sval.npar", "20", target_const_ar)
    println("  ")
    
    push!(modpar.parar, (parmin=0.0 , parmax=10.0, npar=10, nrepeat=5))
    push!(modpar.parar, (parmin=11.0, parmax=21.0, npar=10, nrepeat=5))
    push!(modpar.parar, (parmin=22.0, parmax=32.0, npar=10, nrepeat=5))
    
    return target_const_ar
end #no_pb_x_ginter_no_intra_I0

# This function substitutes the modified variable
function pb_x_ginter_intra(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line("sval.scan_type", "'y'", target_const_ar)
    replace_line("sval.parname",
      "[\"syncoupparI1I2.syn_receptor_par_ar[1].gsyn\",",
      target_const_ar)
    add_line("sval.parname",
      "\"netpar.S_orig[2,1].syn_receptor_par_ar[1].gsyn = " *
      "netpar.S_orig[1,2].syn_receptor_par_ar[1].gsyn\"]",
#     "\"syncoupparI2I1.syn_receptor_par_ar[1].gsyn = " *
#     "syncoupparI1I2.syn_receptor_par_ar[1].gsyn\"]",
      target_const_ar, 1)
    replace_line("sval.parmin", "10.0", target_const_ar)
    replace_line("sval.parmax", "34.0", target_const_ar)
    replace_line("sval.npar", "24", target_const_ar)
    println("  ")
    
    push!(modpar.parar, (parmin=10.0, parmax=20.0, npar=10, nrepeat=5))
    push!(modpar.parar, (parmin=21.0, parmax=31.0, npar=10, nrepeat=5))
    push!(modpar.parar, (parmin=32.0, parmax=35.0, npar=3, nrepeat=5))
    
    return target_const_ar
end #pb_x_ginter_intra

# This function substitutes the modified variable
function no_pb_x_ginter_intra_one_value(original_const_ar, modpar, pardict)

    epsilon = 1.0e-10
    
    target_const_ar = deepcopy(original_const_ar)

    replace_line("sval.parmin", string(pardict["parmin"]), target_const_ar)
    replace_line("sval.parmax", string(pardict["parmax"]), target_const_ar)
    replace_line("sval.npar",   string(pardict["npar"]),   target_const_ar)
    println("  ")
    
    replace_line(r"^syncoupparI1I1.syn_receptor_par_ar\[1\].gsyn",
      string(pardict["gintra"]), target_const_ar)
    replace_line(r"^syncoupparI2I2.syn_receptor_par_ar\[1\].gsyn",
      string(pardict["gintra"]), target_const_ar)

    delg_max_here = pardict["delg_max"] + pardict["gintra"] / 2.0
    
    nfile = ceil(Int64, (1 + delg_max_here - pardict["delg_min"]) /
      pardict["delta_par"] / (pardict["npar_for_one"] + 1))
    println("nfile=", nfile)
	
    for ifile in 1:nfile
        parmin_val = pardict["gintra"] + pardict["delg_min"] + (ifile - 1) *
	  (pardict["npar_for_one"] + 1)* pardict["delta_par"]
	  
	parmax_val_alter = pardict["gintra"] + pardict["delta_par"] *
	  ceil(delg_max_here / pardict["delta_par"])
	
	parmax_val = min(parmin_val + pardict["npar_for_one"] *
	  pardict["delta_par"], parmax_val_alter)
	println("ifile=", ifile, " gintra=", pardict["gintra"], " parmin_val=",
	  parmin_val, " parmax_val=",parmax_val)
	push!(modpar.parar, (parmin=parmin_val, parmax=parmax_val,
	  npar=pardict["npar_for_one"], nrepeat=5))
    end
    
#   push!(modpar.parar, (parmin=20.0, parmax=28.0, npar=8, nrepeat=5))
#   push!(modpar.parar, (parmin=29.0, parmax=37.0, npar=8, nrepeat=5))
#   push!(modpar.parar, (parmin=38.0, parmax=46.0, npar=8, nrepeat=5))
    
    return target_const_ar
end #pb_x_ginter_intra

# This function substitutes the modified variable
function no_pb_x_ginter_intra18(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line("sval.parmin", "20.0", target_const_ar)
    replace_line("sval.parmax", "46.0", target_const_ar)
    replace_line("sval.npar", "26", target_const_ar)
    println("  ")
    
    replace_line(r"^syncoupparI1I1.syn_receptor_par_ar\[1\].gsyn", "18.0",
      target_const_ar)
    replace_line(r"^syncoupparI2I2.syn_receptor_par_ar\[1\].gsyn", "18.0",
      target_const_ar)
      
    push!(modpar.parar, (parmin=20.0, parmax=28.0, npar=8, nrepeat=5))
    push!(modpar.parar, (parmin=29.0, parmax=37.0, npar=8, nrepeat=5))
    push!(modpar.parar, (parmin=38.0, parmax=46.0, npar=8, nrepeat=5))
    
    return target_const_ar
end #pb_x_ginter_intra

# This function substitutes the modified variable
function pb_x_gpb_intra(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line("sval.scan_type", "'y'", target_const_ar)
    replace_line("sval.parname", "[\"cellparI1.gsyn_ext\"]",
      target_const_ar)
    delete_line(
      "\"netpar.S_orig[2,1].syn_receptor_par_ar[1].gsyn = " *
      "netpar.S_orig[1,2].syn_receptor_par_ar[1].gsyn\"]",
#     "\"syncoupparI2I1.syn_receptor_par_ar[1].gsyn = syncoupparI1I2.syn_receptor_par_ar[1].gsyn\"]",
      target_const_ar)
    replace_line("sval.parmin", "0.0", target_const_ar)
    replace_line("sval.parmax", "0.5", target_const_ar)
    replace_line("sval.npar", "20", target_const_ar)
    replace_line(r"^syncoupparI1I2.syn_receptor_par_ar\[1\].gsyn", "16.0",
      target_const_ar)
    replace_line(r"^syncoupparI2I1.syn_receptor_par_ar\[1\].gsyn", "16.0",
      target_const_ar)
    println("  ")
    
    push!(modpar.parar, (parmin=0.0 , parmax=0.2 , npar=10, nrepeat=5))
    push!(modpar.parar, (parmin=0.22, parmax=0.42, npar=10, nrepeat=5))
    push!(modpar.parar, (parmin=0.44, parmax=0.5 , npar=3 , nrepeat=5))
    
    return target_const_ar
end #pb_x_gpb_intra

# This function substitutes the modified variable
function pb_x_gpb_no_inter_no_intra(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line(r"^syncoupparI1I1.syn_receptor_par_ar\[1\].gsyn", "0.0",
      target_const_ar)
    replace_line(r"^syncoupparI1I2.syn_receptor_par_ar\[1\].gsyn", "0.0",
      target_const_ar)
    replace_line(r"^syncoupparI2I1.syn_receptor_par_ar\[1\].gsyn", "0.0",
      target_const_ar)
    replace_line(r"^syncoupparI2I2.syn_receptor_par_ar\[1\].gsyn", "0.0",
      target_const_ar)
    replace_line("netpar.Tper_syn_ext", "200.0", target_const_ar)
    replace_line("netpar.Trand_syn_ext", "10.0", target_const_ar)

    println("  ")
    
    push!(modpar.parar, (parmin=0.0 , parmax=0.2 , npar=10, nrepeat=5))
    push!(modpar.parar, (parmin=0.22, parmax=0.42, npar=10, nrepeat=5))
    push!(modpar.parar, (parmin=0.44, parmax=0.6 , npar=8 , nrepeat=5))
    
    return target_const_ar

end #pb_x_gpb_no_inter_no_intra

# This function substitutes the modified variable
function pb_x_gpb_no_inter_no_intra_s(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line(r"^syncoupparFI1.syn_receptor_par_ar\[1\].gsyn", "1.0",
      target_const_ar)
    replace_line("sval.parmax", "20.0", target_const_ar)
    replace_line("sval.npar", "20", target_const_ar)

    println("  ")
    
    push!(modpar.parar, (parmin=0.0 , parmax=0.2 , npar=10, nrepeat=5))
    push!(modpar.parar, (parmin=0.22, parmax=0.42, npar=10, nrepeat=5))
    push!(modpar.parar, (parmin=0.44, parmax=0.6 , npar=8 , nrepeat=5))
    
    return target_const_ar

end #pb_x_gpb_no_inter_no_intra_s

# This function substitutes the modified variable
function pb_v_inhibition_factor(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line("sval.scan_type", "'y'", target_const_ar)
    replace_line("sval.parname",
      "[\"netpar.inhibition_factor\"]", target_const_ar)
    replace_line("runpar.frac_sd_mnmx", "0.1", target_const_ar)
    
    println("  ")
    
    push!(modpar.parar, (parmin=0.15, parmax=0.55, npar=8, nrepeat=5))
    push!(modpar.parar, (parmin=0.6, parmax=1.0, npar=8, nrepeat=5))

    return target_const_ar
    
end #pb_v_inhibition_factor

# This function substitutes the modified variable
function pb_x_tptbtf_intra(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line("runpar.Tall", "60000.0", target_const_ar)
    replace_line("runpar.tmcol", "70000.0", target_const_ar)
    replace_line("runpar.tstat", "59000.0", target_const_ar)
    replace_line("runpar.traster", "60000.0", target_const_ar)

    println("  ")
    
    return target_const_ar
end #pb_x_tptbtf_intra

# This function substitutes the modified variable
function pb_x_tptbtf_intra_small_gbp(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line("cellparI_mn.gsyn_ext", "0.05", target_const_ar)

    println("  ")
    
    return target_const_ar
end #pb_x_tptbtf_intra_small_gbp

# This function substitutes the modified variable
function pb_tptbtf_x_g_prebot_intra(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line("sval.scan_type", "'y'", target_const_ar)
    replace_line("sval.parmax", "0.5", target_const_ar)
    replace_line("sval.npar", "20", target_const_ar)

    println("  ")
    
    push!(modpar.parar, (parmin=0.0, parmax=0.08, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=0.1, parmax=0.18, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=0.2, parmax=0.28, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=0.3, parmax=0.38, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=0.4, parmax=0.48, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=0.5, parmax=0.5 , npar=0, nrepeat=5))
    
    return target_const_ar
end #pb_tptbtf_x_g_prebot_intra

# This function substitutes the modified variable
function pb_tptbtf_x_Tper_70_intra(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line("sval.parname", "[\"netpar.Tper_syn_ext\"]",
      target_const_ar)
    replace_line("sval.parmin", "100.0", target_const_ar)
    replace_line("sval.parmax", "150.0", target_const_ar)
    replace_line("sval.npar", "20", target_const_ar)
    replace_line("sval.nrepeat", "5", target_const_ar)
    replace_line("netpar.Trand_syn_ext", "10.0", target_const_ar)

    println("  ")
    
    push!(modpar.parar, (parmin=100.0 , parmax=140.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=150.0 , parmax=190.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=200.0 , parmax=240.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=250.0 , parmax=290.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=300.0 , parmax=340.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=350.0 , parmax=390.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=400.0 , parmax=480.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=500.0 , parmax=580.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=600.0 , parmax=640.0, npar=2, nrepeat=5))
    
    return target_const_ar
end #pb_tptbtf_x_Tper_70_intra

# This function substitutes the modified variable
function pb_tptbtf_x_Tper_20_intra(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line("netpar.Tup_syn_ext", "20.0", target_const_ar)

    println("  ")
    
    push!(modpar.parar, (parmin=100.0 , parmax=140.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=150.0 , parmax=190.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=200.0 , parmax=240.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=250.0 , parmax=290.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=300.0 , parmax=340.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=350.0 , parmax=390.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=400.0 , parmax=480.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=500.0 , parmax=580.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=600.0 , parmax=640.0, npar=2, nrepeat=5))

    return target_const_ar
end #pb_tptbtf_x_Tper_20_intra

# This function substitutes the modified variable
function pb_tptbtf_x_Tper_70_intra_Trand_70(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line("netpar.Trand_syn_ext", "70.0", target_const_ar)

    println("  ")
    
    push!(modpar.parar, (parmin=150.0 , parmax=190.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=200.0 , parmax=240.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=250.0 , parmax=290.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=300.0 , parmax=340.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=350.0 , parmax=390.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=400.0 , parmax=480.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=500.0 , parmax=580.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=600.0 , parmax=640.0, npar=2, nrepeat=5))

    return target_const_ar
end #pb_tptbtf_x_Tper_70_intra_Trand_70

# This function substitutes the modified variable
function pb_tptbtf_x_Tper_20_intra_Trand_70(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line("netpar.Trand_syn_ext", "70.0", target_const_ar)

    println("  ")
    
    push!(modpar.parar, (parmin=100.0 , parmax=140.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=150.0 , parmax=190.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=200.0 , parmax=240.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=250.0 , parmax=290.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=300.0 , parmax=340.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=350.0 , parmax=390.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=400.0 , parmax=480.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=500.0 , parmax=580.0, npar=4, nrepeat=5))
    push!(modpar.parar, (parmin=600.0 , parmax=640.0, npar=2, nrepeat=5))

    return target_const_ar
end #pb_tptbtf_x_Tper_20_intra_Trand_70

# This function substitutes the modified variable
function pb_tptbtf_x_Tper_70_intra_Trand_150(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line("netpar.Trand_syn_ext", "150.0", target_const_ar)

    println("  ")
    
    push!(modpar.parar, (parmin=220.0 , parmax=250.0, npar=3, nrepeat=5))
    push!(modpar.parar, (parmin=260.0 , parmax=290.0, npar=3, nrepeat=5))
    push!(modpar.parar, (parmin=300.0 , parmax=330.0, npar=3, nrepeat=5))
    push!(modpar.parar, (parmin=340.0 , parmax=370.0, npar=3, nrepeat=5))
    push!(modpar.parar, (parmin=380.0 , parmax=410.0, npar=3, nrepeat=5))
    push!(modpar.parar, (parmin=420.0 , parmax=450.0, npar=3, nrepeat=5))
    push!(modpar.parar, (parmin=460.0 , parmax=490.0, npar=3, nrepeat=5))
    push!(modpar.parar, (parmin=500.0 , parmax=530.0, npar=3, nrepeat=5))
    push!(modpar.parar, (parmin=540.0 , parmax=570.0, npar=3, nrepeat=5))
    push!(modpar.parar, (parmin=580.0 , parmax=610.0, npar=3, nrepeat=5))
    push!(modpar.parar, (parmin=620.0 , parmax=640.0, npar=2, nrepeat=5))

    return target_const_ar
end #pb_tptbtf_x_Tper_70_intra_Trand_150

# This function substitutes the modified variable
function pb_tptbtf_x_Tper_20_intra_Trand_150(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line("netpar.Trand_syn_ext", "150.0", target_const_ar)

    println("  ")
    
    push!(modpar.parar, (parmin=130.0 , parmax=160.0, npar=3, nrepeat=5))
    push!(modpar.parar, (parmin=170.0 , parmax=200.0, npar=3, nrepeat=5))
    push!(modpar.parar, (parmin=210.0 , parmax=240.0, npar=3, nrepeat=5))
    push!(modpar.parar, (parmin=250.0 , parmax=280.0, npar=3, nrepeat=5))
    push!(modpar.parar, (parmin=290.0 , parmax=320.0, npar=3, nrepeat=5))
    push!(modpar.parar, (parmin=330.0 , parmax=360.0, npar=3, nrepeat=5))
    push!(modpar.parar, (parmin=370.0 , parmax=400.0, npar=3, nrepeat=5))
    push!(modpar.parar, (parmin=410.0 , parmax=440.0, npar=3, nrepeat=5))
    push!(modpar.parar, (parmin=450.0 , parmax=480.0, npar=3, nrepeat=5))
    push!(modpar.parar, (parmin=490.0 , parmax=520.0, npar=3, nrepeat=5))
    push!(modpar.parar, (parmin=530.0 , parmax=560.0, npar=3, nrepeat=5))
    push!(modpar.parar, (parmin=570.0 , parmax=600.0, npar=3, nrepeat=5))
    push!(modpar.parar, (parmin=610.0 , parmax=640.0, npar=3, nrepeat=5))

    return target_const_ar
end #pb_tptbtf_x_Tper_20_intra_Trand_150

# This function substitutes the modified variable
function pb_x_gFIp_no_gPT(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line("sval.parname",
      "[\"netpar.S_orig[3,2].syn_receptor_par_ar[1].gsyn\",",
#     "[\"syncoupparFI2.syn_receptor_par_ar[1].gsyn\"]",
      target_const_ar)
       delete_line("\"netpar.S_orig[2,1].syn_receptor_par_ar[1].gsyn = " *
      "netpar.S_orig[1,2].syn_receptor_par_ar[1].gsyn\"]",
#      delete_line("\"syncoupparI2I1.syn_receptor_par_ar[1].gsyn = syncoupparI1I2.syn_receptor_par_ar[1].gsyn\"]",
       target_const_ar)
    replace_line("sval.parmin", "0.0", target_const_ar)
    replace_line("sval.parmax", "3.0", target_const_ar)
    replace_line("sval.npar", "10", target_const_ar)

    println("  ")
    
    push!(modpar.parar, (parmin=0.0 , parmax=1.6, npar=8, nrepeat=5))
    push!(modpar.parar, (parmin=1.8 , parmax=3.4, npar=8, nrepeat=5))
    
    return target_const_ar
end #pb_x_gFIp_no_gPT

# This function substitutes the modified variable
function pb_x_gFIp_gPT_0_5(original_const_ar, modpar)

    target_const_ar = deepcopy(original_const_ar)

    replace_line("sval.parmin", "0.0", target_const_ar)
    replace_line("sval.parmax", "3.0", target_const_ar)
    replace_line("sval.npar", "10", target_const_ar)
    replace_line("cellparI2.gsyn_ext", "0.5", target_const_ar)
    
    println("  ")
    
    push!(modpar.parar, (parmin=0.0 , parmax=1.6, npar=8, nrepeat=5))
    push!(modpar.parar, (parmin=1.8 , parmax=3.4, npar=8, nrepeat=5))
    
    return target_const_ar
end #pb_x_gFIp_gPT_0_5

# This function generates data files for various levels of gintra.
function generate_files_for_various_gintra(sim_description, array_to_interit, run_dict_str)

#    eval(Meta.parse(compose_write_command(sim_description, array_to_interit,
#      run_dict_str)))
#    println("commamd=", compose_write_command(sim_description,
#      array_to_interit, run_dict_str))

    pardict = Dict()
    pardict["parmin"] = 20.0
    pardict["parmax"] = 46.0
    pardict["npar"] = 26

    pardict["delta_par"] = 1.0
    pardict["delg_min"] = 2.0
    pardict["delg_max"] = 26.0
    pardict["delg_max_slope"] = 0.5
    pardict["npar_for_one"] = 8
    

    no_pb_x_ginter_intra_various(original_const_ar, modpar) =
      no_pb_x_ginter_intra_one_value(original_const_ar, modpar, pardict)
    for xgintra in collect(0.0 : 1.0 : 18.001)
        println("xgintra=", xgintra)
        pardict["gintra"] = xgintra
        write_target_file(no_pb_x_ginter_no_intra_const_ar,
          no_pb_x_ginter_intra_various, "no_pb_x_ginter_intra" *
	  string(xgintra), run_dict_dict["qpd"])
    end

end #generate_files_for_various_gintra

#main

run_here_fl = open("run_here.com", "w")
run_qsub_fl = open("run_qsub.com", "w")
run_qpd_fl = open("run_qpd.com", "w")

num_tasks_here_fl = open("num_tasks.dat", "w")
num_tasks_here_qpd_fl = open("num_tasks_qpd.dat", "w")

here_dict = Dict([("dir", "dir_here"), ("run_fl", run_here_fl),
  ("command", "julia ../prog5/irt.jl dat_here/")])
qsub_dict = Dict([("dir", "dir_qsub"), ("run_fl", run_qsub_fl),
  ("command", "irt_qsub ")])
qpd_dict = Dict([("dir", "dir_qpd"), ("run_fl", run_qpd_fl),
  ("command", "irt_qpd ")])
run_dict_dict = Dict([("here", here_dict), ("qsub", qsub_dict),
  ("qpd", qpd_dict)])

isdir("dir_here") ? nothing : mkdir("dir_here")
isdir("dir_qsub") ? nothing : mkdir("dir_qsub")
isdir("dir_qpd") ? nothing : mkdir("dir_qpd")

# preBot input, V and raster, g_intra. #  c46
dir_orig_file = "pb_v_rast_intra_orig"

#println(runfl, run_command * dir_orig_file)

orig_const_ar = read_org_file(dir_orig_file)

# preBot input, V and raster, g_intra.  c46
eval(Meta.parse(compose_write_command("pb_v_rast_intra",
  "orig_const_ar", "here")))

# preBot input, V and raster, no g_intra.  c461
eval(Meta.parse(compose_write_command("pb_v_rast_no_intra",
  "orig_const_ar", "here")))

# no preBot input, V and raster, g_intra.  c44
eval(Meta.parse(compose_write_command("no_pb_v_rast_intra",
  "orig_const_ar", "here")))
  
# no preBot input, V and raster, no g_intra.  c441
eval(Meta.parse(compose_write_command("no_pb_v_rast_no_intra",
  "no_pb_v_rast_intra_const_ar", "here")))
  
# no preBot input, V and raster, g_intra, asymmetric firing pattern.  c4401
eval(Meta.parse(compose_write_command("no_pb_v_rast_intra_as",
  "no_pb_v_rast_intra_const_ar", "here")))
  
# no preBot input, V and raster, no g_intra, asymmetric firing pattern.  c4411
#  runfl)
eval(Meta.parse(compose_write_command("no_pb_v_rast_no_intra_as",
  "no_pb_v_rast_no_intra_const_ar", "here")))

# no preBot input, x-g_inter, g_intra.  c45
eval(Meta.parse(compose_write_command("no_pb_x_ginter_intra",
  "no_pb_v_rast_intra_const_ar", "qsub")))

# no preBot input, x-g_inter, g_intra.  c452
eval(Meta.parse(compose_write_command("no_pb_x_ginter_intra_I0",
  "no_pb_v_rast_intra_const_ar", "qsub")))
  
# no preBot input, x-g_inter, no g_intra.  c451
eval(Meta.parse(compose_write_command("no_pb_x_ginter_no_intra",
  "no_pb_v_rast_no_intra_const_ar", "qsub")))

# no preBot input, x-g_inter, no g_intra.  c453
eval(Meta.parse(compose_write_command("no_pb_x_ginter_no_intra_I0",
  "no_pb_v_rast_no_intra_const_ar", "qsub")))

# no preBot input, x-g_inter, g_intra6 .
#eval(Meta.parse(compose_write_command("no_pb_x_ginter_intra6",
#  "no_pb_x_ginter_no_intra_const_ar", "qpd")))

# no preBot input, x-g_inter, g_intra18 .
#eval(Meta.parse(compose_write_command("no_pb_x_ginter_intra18",
#  "no_pb_x_ginter_no_intra_const_ar", "qpd")))

# no preBot input, x-g_inter, various levels of g_intra.
generate_files_for_various_gintra("no_pb_x_ginter_intra_various",
  "no_pb_x_ginter_no_intra_const_ar", "qpd")

# preBot input, x-g_inter, g_intra.  c47
eval(Meta.parse(compose_write_command("pb_x_ginter_intra",
  "orig_const_ar", "qsub")))

# preBot input, x-g_inter, g_intra.  c471
eval(Meta.parse(compose_write_command("pb_x_gpb_intra",
  "pb_x_ginter_intra_const_ar", "qsub")))

# preBot input, x-g_inter, no g_inter, no g_intra.
eval(Meta.parse(compose_write_command("pb_x_gpb_no_inter_no_intra",
  "pb_x_gpb_intra_const_ar", "qsub")))

# preBot input, x-g_inter, no g_inter, no g_intra, g_Fr=1mS/cm^2.
eval(Meta.parse(compose_write_command("pb_x_gpb_no_inter_no_intra_s",
  "pb_x_gpb_no_inter_no_intra_const_ar", "qsub")))

# preBot input, x_inhibition_factor.
eval(Meta.parse(compose_write_command("pb_v_inhibition_factor",
  "orig_const_ar", "qsub")))

# preBot input, tptbtf, g_intra.  c43
eval(Meta.parse(compose_write_command("pb_x_tptbtf_intra",
  "orig_const_ar", "here")))

# preBot input, tptbtf, g_intra.  c431
eval(Meta.parse(compose_write_command("pb_x_tptbtf_intra_small_gbp",
  "pb_x_tptbtf_intra_const_ar", "here")))

# preBot input, tptbtf, x-g_prebot, g_intra.  c4341
eval(Meta.parse(compose_write_command("pb_tptbtf_x_g_prebot_intra",
  "pb_x_tptbtf_intra_const_ar", "qsub")))

# preBot input, tptbtf, Tup=70, x-g_prebot, g_intra.  c43301
eval(Meta.parse(compose_write_command("pb_tptbtf_x_Tper_70_intra",
  "pb_tptbtf_x_g_prebot_intra_const_ar", "qsub")))

# preBot input, tptbtf, Tup=20, x-g_prebot, g_intra.  c43201
eval(Meta.parse(compose_write_command("pb_tptbtf_x_Tper_20_intra",
  "pb_tptbtf_x_Tper_70_intra_const_ar", "qsub")))

# preBot input, tptbtf, Tup=70, x-g_prebot, g_intra.
eval(Meta.parse(compose_write_command("pb_tptbtf_x_Tper_70_intra_Trand_70",
  "pb_tptbtf_x_Tper_70_intra_const_ar", "qsub")))

# preBot input, tptbtf, Tup=20, x-g_prebot, g_intra.
eval(Meta.parse(compose_write_command("pb_tptbtf_x_Tper_20_intra_Trand_70",
  "pb_tptbtf_x_Tper_20_intra_const_ar", "qsub")))

# preBot input, tptbtf, Tup=70, x-g_prebot, g_intra.
eval(Meta.parse(compose_write_command("pb_tptbtf_x_Tper_70_intra_Trand_150",
  "pb_tptbtf_x_Tper_70_intra_const_ar", "qsub")))

# preBot input, tptbtf, Tup=20, x-g_prebot, g_intra.
eval(Meta.parse(compose_write_command("pb_tptbtf_x_Tper_20_intra_Trand_150",
  "pb_tptbtf_x_Tper_20_intra_const_ar", "qsub")))

# preBot input,x-gFIp, gPT=0, g_intra.  c522
eval(Meta.parse(compose_write_command("pb_x_gFIp_no_gPT",
  "pb_x_ginter_intra_const_ar", "qsub")))

# preBot input,x-gFIp, gPT=0.5, g_intra.  c511
eval(Meta.parse(compose_write_command("pb_x_gFIp_gPT_0_5",
  "pb_x_gFIp_no_gPT_const_ar", "qsub")))

close(run_here_fl)
close(run_qsub_fl)
close(run_qpd_fl)
close(num_tasks_here_fl)
close(num_tasks_here_qpd_fl)
