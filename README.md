# out-of-equilibrium-interlayer-exchange-coupling-
Evaluates the out-of-equilibrium interlayer exchange coupling for multi-layer systems with either a single or resonant tunnel barrier in a 2-band, tight-binding model.

Systems of the form LEAD/INS/FM/NM/FM/LEAD are considered. LEAD's are conducting semi-infinite, FM's are ferromagnetic materials, NM is a conducting non-magnetic spacer and INS is either a single insulating layer (ooeIEC_single.cpp) or a resonant tunnel barrier of the form INS/NM/INS (ooeIEC_double.cpp). 

Parameters specifically model Cu/INS/Co/Cu/Co/Cu where the insulating layers are modelled using the minority spin Co parameters.

Layer thicknesses (atomic planes) are specified globally except FM tr-layer spacer which is considered at thicknesses 2,4,6,8,10.

Bias, temperature, fermi-level, imaginary part of energy and several other system parameters also declared globally.

Primary function "f" evaluates the surface Greens functions up to a clevage plane set in the center of the NM between FM's and returns spin current density functions (to be integrated) for either the IEC or STT terms depending on the value of global flag. Theta integration of STT term also carried out in "f".

"eqm_func" evalutes the final IEC term integral (itself has two seperate terms for left and right Fermi levels) using basic weighted sum for k-space and sum over energy residues (Matsubara). Energy integral optimised enough to allow for k-space integral to be evaluted multiple times at increasing grid density to check for convergence.

"oee_func" evalutes the final STT term integral using fixed grid weighted sum in k-space and NAG adaptive integration routine for energy. 

Could be easily converted for different material system configuration (provided leads semi-infinite) and larger number of bands.
