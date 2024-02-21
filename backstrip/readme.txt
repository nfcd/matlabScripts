Explanation:

Backstrip performs 1D backstripping of sedimentary strata.

Assumptions are: Exponential decrease of porosity with depth
                 Fully Saturated sediments
                 Local compensation of sedimentary loads (zero flexural strength)
                 

Program structure:

The program is designed such that the user follows a logical sequence of steps. I have
intentionally introduced some "flags" that make sure the user follows that sequence. The path
of solution is:

PARAMETERS

Define tops and bases of units (present geometry)
Define ages
Define densities of water and mantle, and dry densities of sediments
Define porosities and porosity depth coefficients of sediments

The user can input these parameters through a window , or she can load these 
parameters from a txt (tab delimited) file with a number of lines equal to the number
of sedimentary units, each line having:


top 
base 
age of_top 
age of_base 
dry density 
por depth coeff 
surface porosity 
basin setting (1 if the basin is marine. 2 if the basin is continental)

For an example see file insample.txt

SOLVE

1. Decompact sediments
2. Correct for the weight of the sediments
3. Correct for sea level changes and water depth: This algorithm has improved thanks 
   to Elias Gomez, who pointed out that relative instead of absolute water changes must 
   be used in the corrections. SEA LEVELS AND WATER DEPTHS CAN ONLY BE ENTERED THROUGH A
   DIALOGUE WINDOW.


PLOTS

1. Decompacted units vs. time
2. Total thickness, total decompacted thickness, total decompacted thickness corrected
   for the weight of sediments   vs. time
3. Tectonic subsidence


Plots can be exported or saved using the matlab menu

RUNNING BACKSTRIP

In the matlab command window simply change to the backstrip folder

Then type backstrip

ANY QUESTIONS EMAIL nfcd@mac.com

These scripts are copyright of Nestor Cardozo 2009. They can only be used for academic/research purposes and they cannot be distributed by third parties

Nestor Cardozo assumes no liability for damages, direct or consequential, which may result from the use of the scripts