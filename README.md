# multistrainSIR
Simulation Code for "The local stability of a modified multi-strain SIR model for emerging viral strains" by Fudolig and Howard submitted to PLOS One.

simulation.R simulates the code for reproducing the images in the manuscript. simulation.R needs the following packages:
-ggplot2
-tidyR
-gridExtra

The code also defines two functions:

-for_plots produces the surveillance plots given the transmission and removal coefficients for both strains
-surv produces the colormap of the equilibrium points reached for a range of values of the reproductive numbers for both strains.
