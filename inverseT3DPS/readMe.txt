inverse: A Matlab group of scripts to find the true 3D trishear model that best restores a folded bed surface to a plane

This optimized 3D trishear inversion uses a pattern search constrained optimization Matlab function (pattern search). The inversion is much faster than the standard grid search inversion, and is not affected as much by local minima.

true 3D trishear algorithm from Cristallini et al (2004): GSA bulletin 116, 938-952

Application: 

In Matlab change directory to the inverseT3D folder

To illustrate the use of these scripts, an example from a synthetic true 3D trishear model is included. This model has a ramp angle of 30 degrees, a P/S varying linearly from 1.0 to 2.0 along the tip line, a trishear angle varying linearly from 40.0 to 80.0 degrees along the tip line, slip of 100 units, and slip rake of 90 degrees (slip is perpendicular to fault strike)

The geometry of one of the beds of the true 3D synthetic model is included in file bed.txt. This file contains the x, y, and z coordinates of points on the bed surface

File invBed.m shows the application of the scripts to the included bed. For this case we will assume than the P/S, trishear angle, and slip at both fault tips are unknown. The search will be then for six parameters. Type invBed and watch the program search for the parameters that best restore the bed to a plane using a Matlab constrained optimization function (fmincon). When the search is finished, the program plots the point data and the best fit surface

To run an inversion you need:

1. A text file containing the x, y and z position of points on the bed surface to be restored

2. An initial estimate of the parameters to be searched. CAREFUL: This estimate has to be within the limits of your search

How does the program work?

The file invBed.m shows the application of the scripts to the included bed. This file reads the points on the bed, sets the initial estimate of the parameters (initial guess), runs the inversion, and plots the results. In a file like this is where you set the initial estimate of the parameters (the initial guess). CAREFUL: The initial estimate has to be within the limits of your search

Function gradcongen.m runs the inversion and gives the results in a structure called history. After the inversion, type history.x to see the history of parameter values, or history.fval to see the history of the objective function through the search. The last line of history.x corresponds to the best fit parameters, that combination of parameters for which the objective function is minimum. In gradcongen.m is where you decide which parameters to search and the limits of these parameters. CAREFUL: Make sure that the searched parameters are of the same order of magnitude. For that you might need to scale the parameters before the search

Function backtrishear.m computes the error between the restored bed and the plane that best fits the restored bed. This is the objective function. The aim of the inversion is to find the combination of trishear parameters that minimizes the objective function

The following two functions are used for plotting the best fit model:

restorebed.m : Restores the bed and finds the the plane that best fits the restored bed (the best fit plane)

deformbed.m : Deforms the best fit plane. This deformed best fit plane is the best fit trishear model. This is usually plotted over the data to visually check how good the fit is

Function veltrishear.m is the velocity field of the true 3D trishear model as described in Cristallini et al. (2004). For a discussion of the problems related to this formulation check Cardozo (2008): Journal of Structural Geology 30, 327-340

NOTE: Slip should be positive for reverse faults, and negative for normal faults. The P by S is positive in both reverse and normal faults

These scripts are copyright of Nestor Cardozo 2009. They can only be used for academic/research purposes and they cannot be distributed by third parties

Nestor Cardozo assumes no liability for damages, direct or consequential, which may result from the use of the scripts