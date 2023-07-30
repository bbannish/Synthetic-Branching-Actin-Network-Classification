# Growth of 2D synthetic actin network

Summary

We provide code associated with the growth of synthetic dendritic actin network. The two-dimensional model is developed in Matlab for a general cell type but incorporating specific protein and protein interactions is straightforward. The model consists of discrete stochastic processes associated with the assembly/disassembly, capping, and branching of actin filaments with biologically informed reaction rates. 

To run the model

To run “stochastic_code_serial.m”, simply open it in Matlab and click the “Run” button. The code is currently set up to do 100 independent simulations of the “Label 5” conditions: unlimited Arp2/3, limited monomers, and capping with capping probability 0.00025. To change the number of independent simulations, edit the “run_number = 100” line to be whatever number of simulations is desired. By setting various combinations of “capyes”, “Arp23yes”, and “monomeryes” equal to 1 (include that feature in the run) or 0 (do not include that feature in the run), the other growth conditions can be achieved. The 0.05 constant that is currently in the definition of “pcap” can be changed to get the various capping probabilities. The result of running the code (as is) is a movie of network growth saved as “network_movie.avi”, a figure of the final network (which corresponds to Figure 2g from the paper), a gray-scale image such as those from Figure 4, and “Factin_array_cap_monomer.mat” (a matrix of the actin densities at the final time point, T=10. Each column corresponds to an independent simulation. The 49x49 density data from a single run is reformatted into a 2401x1 vector to more easily fit in “Factin_array_cap_monomer.mat”.) 

Paper

Inferring local molecular dynamics from the global actin network structure: a case study of 2D synthetic branching actin networks

Minghao W. Rostami, Brittany E. Bannish, Kelsey Gasior, Rebecca L. Pinals, Calina Copos and Adriana T. Dawes

bioRxiv: https://www.biorxiv.org/content/10.1101/2022.09.06.506753v2

doi: https://doi.org/10.1101/2022.09.06.506753
