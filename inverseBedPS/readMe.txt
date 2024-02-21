inverseBed: A Matlab group of scripts to find the combination of trishear parameters that best restores a folded bed to a straight line

This is the classical trishear inversion as described in Allmendinger (1998): Tectonics 17, 640-656. The only difference is that the inversion is optimized using a pattern search constrained optimization Matlab function (patternsearch). The inversion is therefore much faster than the standard grid search inversion, and is not affected as much by local minima.

Application: 

In Matlab change directory to the inverseBed folder

To illustrate the use of these scripts, an example from a trishear like fault propagation fold in West Central Taiwan (Lin et al. 2007. Journal of Structural Geology 29, 1267-1280) is included

The geometry of the beds in this structure is defined by text files, one for each bed, and each one containing the x y coordinates of points on the bed. See files bed5.txt, bed6.txt, bed7.txt, bedGrowth4.txt, bedGrowth3.txt, bedGrowth2.txt, bedGrowth1.txt, and ground.txt. To plot the beds, type plotBeds in the Matlab command window. You will see the pregrowth beds in black and the growth beds in blue, and a red rectangle which is the area where the final location of the fault tip will be searched

Files invBed5.m, invBed6.m, and invBed7.m show the application of the scripts to beds 5, 6,and 7, respectively. Type invBed5.m and watch the program search for the parameters that best fit bed 5. When the search is finished, the program plots the data in red and the best fit model in blue. Run files invBed6.m and invBed7.m as well

Once you understand the scripts, do the following exercise: 

Use the tip location, PS, trishear angle, and ramp angle from the inversion of the pregrowth strata 5, 6, and 7, and assume these parameters to be fixed. Then, search for the slip that best restores the growth strata. You will restore progressively the structure using the growth strata as Lin et al. (2007) did. Compare your results with these authors.

To run an inversion you need:

1. A text file containing the x y position of points on the bed to be restored

2. An initial estimate of the parameters to be searched. CAREFUL: This estimate has to be within the limits of your search

How does the program work?

The file invBed5.m, invBed6.m, or invBed7.m shows the application of all scripts for the example of Lin et al. (2007). This file reads the points on the bed, sets the initial estimate of the parameters (initial guess), runs the inversion, and plots the results. In a file like this is where you set the initial estimate of the parameters (the initial guess). CAREFUL: The initial estimate has to be within the limits of your search

Function gradcongen.m runs the inversion and gives the results in a structure called history. After the inversion, type history.x to see the history of parameter values, or history.fval to see the history of the objective function through the search. The last line of history.x corresponds to the best fit parameters, that combination of parameters for which the objective function is minimum. In gradcongen.m is where you decide which parameters to search and the limits of these parameters. CAREFUL: Make sure that the searched parameters are of the same order of magnitude. For that you might need to scale the parameters before the search

Function backtrishear.m computes the error between the restored bed and the straight line that best fits the restored bed. This is the objective function. The aim of the inversion is to find the combination of trishear parameters that minimizes the objective function. For a better explanation of the inversion strategy, read Richard Allmendinger (1998)

The following two functions are used for plotting the best fit model:

restorebed.m : Restores the bed and finds the line that best fits the restored bed (the best fit line)

deformbed.m : Deforms the best fit line. This deformed best fit line is the best fit trishear model. This is usually plotted over the actual data to visually check how good the fit is

Function veltrishear.m is the velocity field of the trishear model. The program uses the simplest velocity field of trishear: A linear in x velocity field (Zehnder and Allmendinger, 2000) 

NOTE: Slip and slip increment should be positive for reverse faults, and negative for normal faults. The P by S is positive in both reverse and normal faults

These scripts are copyright of Nestor Cardozo 2009. They can only be used for academic/research purposes and they cannot be distributed by third parties

Nestor Cardozo assumes no liability for damages, direct or consequential, which may result from the use of the scripts