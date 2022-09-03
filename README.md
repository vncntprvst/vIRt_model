# vIRt
Simulation programs for the article of Golomb et al., 2022.

This github site includes computer programs and scripts for generating figures
from the article:

David Golomb, Jeffrey D. Moore, Arash Fassihi, Jun Takatoh, Vincent Prevosto,
Fan Wang and David Kleinfeld,
Theory of hierarchically-organized neuronal oscillator dynamics that mediate
rodent rhythmic whisking.
Neuron, in press, 2022.

Software needed:
julia compiler, 
XMGrace
LaTeX

# Generate data files used for the figures in the article:

cd datfig

julia ../genfig/gen_input_files.jl

# Generate figure 6A,B (and more):

julia ../../prog/irt.jl pb_v_rast_intra

julia ../../prog/irt.jl pb_v_rast_no_intra

../genfig/scripts_fig/vtfig_b.com

../genfig/scripts_fig/rasit_b.com

# Generate figure S2:

[Sumulation of the conductance-based model]

cd datfig

julia ../prog/irt.jl no_pb_x_ginter_no_intra

julia ../prog/irt.jl no_pb_x_ginter_no_intra_I0

[Simulation of the rate model]

cd ../rate_num/datfig

julia ../prog_rf/rf.jl no_pb_x_ginter_no_intra

julia ../prog_rf/rf.jl no_pb_x_ginter_no_intra_I0

[Analytical solution]

julia ../../rate_ode/tvnc.jl no_pb_x_ginter_no_intra

julia ../../rate_ode/tvnc.jl no_pb_x_ginter_no_intra_I0

cd ../../datfig/

[generating the figure]

../genfig/scripts_fig/cmta_c.jl

[One needs to save the eps files from XMGrace]

cd ../figs2_gen/

latex figS2.tex

dvipdf figS2

evince figS2.pdf

