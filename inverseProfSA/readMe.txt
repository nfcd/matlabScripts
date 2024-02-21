inverseProf: A Matlab program to fit a 2D trishear model to bedding intersections along a topographic profile

This is the inversion described in Cardozo (2005), Journal of Structural Geology 27, 495-502. The only difference is that the inversion is optimized using a simulated annealing, constrained optimization Matlab function (simulannealbnd). The inversion is therefore much faster than the standard grid search inversion, and is not affected as much by local minima.

Use: 

In Matlab change directory to the inverseProf folder

The file invprof.m is an example of the use of the program for the situation shown in figure model.jpg. We assume the synthetic to be the "unknown" model, and the bedding intersections of the synthetic model with the topographic profile (blue line in model.jpg) the data. We search for the combination of all six trishear parameters (x,y tip location, ramp angle, p by s, trishear angle, and slip) that best reproduces the bedding intersections along the topographic profile. If the inversion works, the best fit model should be similar to model.jpg 

Type invprof in the Matlab command window. You will see how the program searches for the combination of trishear parameters that best fit the bedding intersections along the topographic profile of model.jpg. When finished, the program plots the profile, the bedding intersections, and the best trishear model

Input files:

To run an inversion you need:

1. A text file with the x and y positions of the points on the topographic profile. For an example see file profile.txt (topographic profile of model.jpg) 

2. A text file with the x, y, and dip (degrees) of the intersections. Positive dips are towards the left and negative dips towards the right. For an example see file intersections.txt (bedding intersections along topographic profile of model.jpg) 

3. A text file with the stratigraphy in an area outside the fold in the hanging wall or footwall. This stratigraphy should correspond with the stratigraphy of the contacts along the profile. This file should contain the following lines: 

- A first line corresponding to the x location where the stratigraphy is measured. 

- n lines (n being the number of contacts or layers) with the y position (elevation) of the contacts, starting from the oldest contact (lowest y) and increasing upwards

- A line with the regional dip of the contacts. Positive dip is towards the left and negative dip towards the right

- A line with an number that can be 0 or 1, to indicate where the stratigraphy is measured: 0 if in the footwall, or 1 if in the hanging wall

For an example see file undstrat.txt (stratigraphy in hanging wall area outside the fold of model.jpg)

4. An initial estimate of the parameters to be searched. CAREFUL: This estimate has to be within the limits of your search

How does the program work?

The file invprof.m is where you plan the inversion: Read the input data, set the initial estimate, run the inversion, and plot the results

The function gradcongen.m runs the inversion and gives the results in a structure called history. After the inversion, type history.x to see the history of parameter values, or history.fval to see the history of the objective function through the search. The last line of history.x corresponds to the best fit parameters, that combination of parameters for which the objective function is minimum. In gradcongen.m is where you define the parameters to be searched and their limits. CAREFUL: Make sure that the searched parameters are in the same order of magnitude. For that you might need to scale the parameters before the search

The function backtrishear.m computes the error between the actual and the modeled intersections. This is the objective function. The aim of the inversion is to find the combination of trishear parameters that minimize the objective function. Modify this file if you need to compute the error (objective function) in a different way

The function deformbed.m runs the trishear model forward

The function veltrishear.m is the velocity field of the trishear model. The program uses the simplest velocity field of trishear: A linear in x velocity field (Zehnder and Allmendinger, 2000)

NOTE: Slip and slip increment should be positive for reverse faults, and negative for normal faults. The P by S is positive in both reverse and normal faults

These scripts are copyright of Nestor Cardozo 2009. They can only be used for academic/research purposes and they cannot be distributed by third parties

Nestor Cardozo assumes no liability for damages, direct or consequential, which may result from the use of the scripts