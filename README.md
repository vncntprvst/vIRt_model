# vIRt
Simulation programs for the article of [Golomb et al., 2022](https://doi.org/10.1016/j.neuron.2022.08.020).

This github site includes computer programs and scripts for generating figures
from the article:
```
David Golomb, Jeffrey D. Moore, Arash Fassihi, Jun Takatoh, Vincent Prevosto,
Fan Wang and David Kleinfeld,
**Theory of hierarchically-organized neuronal oscillator dynamics that mediate
rodent rhythmic whisking.**
*Neuron*, 110(22):3833-3851.e22, 2022. 
DOI: https://doi.org/10.1016/j.neuron.2022.08.020
```

Software needed:
+ julia compiler 
+ XMGrace
+ LaTeX

### Install Julia environment and packages

1. **Clone the Repository:**
```sh
git clone https://github.com/david-golomb/vIRt.git
```
2. **Activate the project**
Open a Julia terminal and navigate to the project directory:
`cd vIRt` then `julia`.
```julia
import Pkg
Pkg.activate(@__DIR__)
```
3. **Install Required Packages**
```julia
Pkg.add("Formatting")
Pkg.add("SmoothingSplines")
Pkg.add("Polynomials")
Pkg.add("DSP")
Pkg.add("NLsolve")
Pkg.add("StatsBase")
Pkg.add("DataStructures") 
```
4. **Ensure all dependencies are installed**
Now that the Project.toml file exists, run the following commands again:
```julia
Pkg.activate(@__DIR__)
Pkg.instantiate()
```
5. **Test the installation**
Load the Formatting package
```julia
using Formatting
```
then exit the Julia terminal: `exit()`.

### Generate data files used for the figures in the article:
```sh
cd datfig
julia ../genfig/gen_input_files.jl
```
### Generate figure 6A,B (and more):
```sh
cd datfig
julia --project=.. ../prog/irt.jl pb_v_rast_intra_orig
```
```sh
cd ..\datfig\dir_here\
julia --project=../.. ../../prog/irt.jl pb_v_rast_no_intra
```
```sh
cd datfig
../genfig/scripts_fig/vtfig_b.com
```
```sh
cd datfig
../genfig/scripts_fig/rasit_b.com
```

For the two commands above, if on Windows, you can use Git Bash or any other terminal that supports Unix-like commands. Here is how you can do it:

* Open Git Bash or another terminal that supports shell scripts.
* Navigate to the directory containing your script.
* Run the scripts using the sh command.
```sh
sh ../genfig/scripts_fig/vtfig_b.com
sh ../genfig/scripts_fig/rasit_b.com
```

### Generate figure S2:

[Simulation of the conductance-based model]
```sh
cd datfig
julia --project=.. ../prog/irt.jl no_pb_x_ginter_no_intra
julia --project=.. ../prog/irt.jl no_pb_x_ginter_no_intra_I0
```
[Simulation of the rate model]
```sh
cd ../rate_num/datfig
julia ../prog_rf/rf.jl no_pb_x_ginter_no_intra
julia ../prog_rf/rf.jl no_pb_x_ginter_no_intra_I0
```
[Analytical solution]
```sh
julia ../../rate_ode/tvnc.jl no_pb_x_ginter_no_intra
julia ../../rate_ode/tvnc.jl no_pb_x_ginter_no_intra_I0
cd ../../datfig/
```
[generating the figure]
```sh
../genfig/scripts_fig/cmta_c.jl
```
[One needs to save the eps files from XMGrace]
```sh
cd ../figs2_gen/
latex figS2.tex
dvipdf figS2
evince figS2.pdf
```
