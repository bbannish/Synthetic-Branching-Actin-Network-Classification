# actin

README for “stochastic_code_serial.m”

To run “stochastic_code_serial.m”, simply open it in Matlab and click the “Run” button.

The code is currently set up to do 100 independent simulations of the “Label 5” conditions: unlimited Arp2/3, limited monomers, and capping with capping probability 0.00025. 

To change the number of independent simulations, edit the “run_number = 100” line to be whatever number of simulations is desired.

By setting various combinations of “capyes”, “Arp23yes”, and “monomeryes” equal to 1 (include that feature in the run) or 0 (do not include that feature in the run), the other growth conditions can be achieved. The 0.05 constant that is currently in the definition of “pcap” can be changed to get the various capping probabilities.

The result of running the code (as is) is a movie of network growth saved as “network_movie.avi”, a figure of the final network (which corresponds to Figure 2g from the paper), a gray-scale image such as those from Figure 4, and “Factin_array_cap_monomer.mat” (a matrix of the actin densities at the final time point, T=10. Each column corresponds to an independent simulation. The 49x49 density data from a single run is reformatted into a 2401x1 vector to more easily fit in “Factin_array_cap_monomer.mat”.) 
