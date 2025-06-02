# PartIIIProject
 This folder contains but is not limited to the data referenced in the report, and the code used.


 In experiments:
- Ergodic - code to establish l-spin moves ergodic to the single spin flip.
- singlefliptau - code to generate single flip energy autocorrelation times
- singleflip - code to generate energy anneals of single flip system, and compare them to analytical results.

 In Simulations:
-  compare_ac - Generating autocorrelation times for spin, magnetisation and energy for single flip.
-  energy - generating energy from l-spin anneals for a wide range of l
- energy_for_saddles - generating l-spin anneals for ergodic pairs of (L,l) to compare with topological crossover data
- energy_l_p - comparing how l-spin anneals change with varying tau
- l_autocorrelation - generating unnormalised energy, magnetisation and spin autocorrelation functions for T > Tc. normalising, averaging and extracting autocorrelation times.
- saddles - generating saddles data, processing to get histograms, sigmoid fits etc
- saddles_l=1 - generating saddles data just for l = 1
- saddles_scaling - generating saddles data for wide range of L, combining with data from "saddles" to perform finite-size analysis 
- sf_autocorrelation - code to generate single flip energy autocorrelation times, this time saving the autocorrelation functions as csvs on the way.
- VFT_close_to_Tc - attempted VFT and CSD fits from l-spin spin autocorrelation data generated close to Tc. Inconclusive results
- VFT_spin - attempted VFT and CSD fits for the l-spin spin autocorrelation data
- VFT_spin - attempted VFT and CSD fits for the l-spin energy autocorrelation data



