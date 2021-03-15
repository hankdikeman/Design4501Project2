# Design4501Project2
Description: This is a reactor network simulation for the production of OME3-5 from renewable feedstocks (CO2, H2) using the aqueous formaldehyde route. This simulation was performed for the preliminary process evaluation required for report 2 for Design 4501W at the University of Minnesota-Twin Cities. The program evaluates each process stream before and after calculation using material balance models, and adjusts process streams for subsequent iterations using the finite-difference gradient. This process is repeated until the convergence criterion is reached.

How to run: To run this program, install the required dependencies if needed and run the process configuration driver file using the console command 'python3 NetworkDriver.py'

Contributors:
Henry Dikeman
Patrick Gibbons-Peterson
Sandhya Appiah
Shivam Yadav

Program Dependencies:
Python:     3.8.5<=
Numpy:      1.20.1<=
Sympy:      1.7.1<=
Scipy:      1.6.1<=

note: it is likely that backwards compatability with older package versions exists, but these are the versions known to work with the current version of the simulation
