# This file processes the rastergram.

if length(ARGS) >= 1
    suffix = ARGS[1]
else
    suffix = "a1"
end

if length(ARGS) >= 2
    tshift = parse(Float64, ARGS[2])
else
    tshift = 0.0
end
println("tshift=", tshift)

if length(ARGS) >= 3
    ipar = ARGS[3]
else
    ipar = ""
end

if length(ARGS) >= 4
    irepeat = ARGS[4]
else
    irepeat = ""
end

dir_prog = "../../prog/"
dir_dat = pwd() * "/" # "../dat5/"
include(dir_prog * "irt_types.jl")
include(dir_dat * suffix * "/irt_const.jl")

#pop_ch_ar = ["I", "P", "F"]
#pop_ch_ar = ["I", "F"]
pop_ch_ar = netpar.ar_cell_str
Npop = length(pop_ch_ar)

frst_name = suffix * "/irt.rst"
if ipar != ""
   frst_name = frst_name * "." * ipar
   if irepeat != ""
        frst_name = frst_name * "." * irepeat 
   end	
end

frst = open(frst_name, "r")
#frsl = Array{IOStream}(Npop)
frsl = Array{IOStream,1}(undef, Npop)
for ipop in (1: Npop)
    fname = suffix * "/process." * pop_ch_ar[ipop]
    println("fname=", fname)
    frsl[ipop] = open(fname, "w")
end

for line in eachline(frst)
    line_split = split(line, " ")
    ipop = parse(Int64, line_split[3])
    println(frsl[ipop], parse(Float64, line_split[1]) - tshift, " ",
      line_split[4])
end


close(frst)
for ipop in (1: length(pop_ch_ar))
    close(frsl[ipop])
end
