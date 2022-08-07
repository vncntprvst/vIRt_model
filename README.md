# vIRt
Simulation programs for the article of Golomb et al., 2022.

This github site includes computer programs and scripts for generating figures
from the article:

David Golomb, Jeffrey D. Moore, Arash Fassihi, Jun Takatoh, Vincent Prevosto,
Fan Wang and David Kleinfeld,
Theory of hierarchically-organized neuronal oscillator dynamics that mediate
rodent rhythmic whisking.
Neuron, accepted for publication, 2022.

Software needed:
julia compiler
XMGrace

Generate data files used for the figures in the article:
cd datfig
julia ../genfig/gen_input_files.jl

Generate figure 6A,B (and more):

julia ../../prog/irt.jl pb_v_rast_intra

julia ../../prog/irt.jl pb_v_rast_no_intra

../genfig/scripts_fig/vtfig_b.com

../genfig/scripts_fig/rasit_b.com
