2D inverse trishear fault propagation folding with a circular fault and inclined shear in the backlimb. The maximum central angle of the fault is entered by the user. Beyond that, the fault propagates at a constant dip angle.

This is the classical trishear inversion as described in Allmendinger (1998): Tectonics 17, 640-656. The only difference is that the inversion is optimized using a simulated annealing  Matlab function (simulannealbnd in the Matlab Global Optimization Toolbox). The inversion is therefore much faster than the standard grid search inversion, and is not affected as much by local minima.

Application: 

In Matlab change directory to the inverse folder

To illustrate the use of the scripts, the uppermost four beds of the default forward model (trishear.m) in the forward folder are included

The geometry of the beds is defined by text files, one for each bed, and each one containing the x y coordinates of points on the bed. See files bed4.txt, bed5.txt, bed6.txt, and bed7.txt. To plot the beds, type plotBeds in the Matlab command window. You will see the beds in black and a red rectangle which is the area where the center of curvature of the fault is searched

Files invBed4.m, invBed5.m, invBed6.m, and invBed7.m show the application of the scripts to beds 4, 5, 6, and 7, respectively. Type invBed4.m and watch the program search for the parameters that best fit bed 3. When the search is finished, the program plots the data in red and the best fit model in blue. Type history.x to see the history of parameter values, or history.fval to see the history of the objective function through the search. The last line of history.x corresponds to the best fit parameters, that combination of parameters for which the objective function is minimum. 

Try running files invBed5.m, invBed6.m, and invBed7.m as well 

To run an inversion you need:

1. A text file containing the x y position of points on the bed to be restored

2. An initial estimate of the parameters to be searched. CAREFUL: This estimate has to be within the limits of your search

How does the program work?

The files invBed4.m, invBed5.m, invBed6.m or invBed7.m show the application of the scripts for a synthetic model. This file reads the points on the bed, sets the initial estimate of the parameters (initial guess), runs the inversion, and plots the results. In a file like this is where you set the initial estimate of the parameters (the initial guess).

Function gradcongen.m runs the inversion and gives the results in a structure called history.  In gradcongen.m is where you decide which parameters to search and the limits of these parameters. CAREFUL: Make sure that the searched parameters are of the same order of magnitude. For that you might need to scale the parameters before the search

Function backtrishear.m computes the error between the restored bed and a straight line that best fits the restored bed. This is the objective function. The aim of the inversion is to find the combination of trishear parameters that minimizes the objective function. For a better explanation of the inversion strategy, read Allmendinger (1998)

The following two functions are used for plotting the best fit model:

restorebed.m : Restores the bed and finds the straight line that best fits the restored bed (the best fit line)

deformbed.m : Deforms the best fit line. This deformed best fit line is the best fit trishear model. This is usually plotted over the actual data to check how good the fit is

Function veltrishear.m is the velocity field of the trishear model. The program uses the simplest velocity field of trishear: A linear in x velocity field (Zehnder and Allmendinger, 2000) 

NOTE: The shear angle is measured from the vertical. Positive shear angle corresponds to antithetic shear. Slip and slip increment should be positive for reverse faults, and negative for normal faults. The P by S is positive in both reverse and normal faults

These scripts are copyright of Nestor Cardozo 2012. Nestor Cardozo assumes no liability for damages, direct or consequential, which may result from the use of the scripts