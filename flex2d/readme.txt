Flex2d: Computes the deflection of a continuous or broken elastic lithosphere
underlying a fluid astenosphere (analogy: a beam on an elastic foundation). 

The solution algorithm is from Hetenyi (1946).

The program is designed such that the user follows a logical sequence of steps. I have
intentionally introduced "flags or warnings" to follow this sequence. 

The solution steps are:

1. Define elastic parameters of lithosphere

2. Define type of analysis: Continuous or broken lithosphere. Define density of the
   astenosphere

3. Define main geometry: Extent of the region with loads, extent of the region under
   analysis, maximum height of the loads, and number of intervals used to discretize
   the load. For the example in file insample.txt the extent of the loaded region is
   240 km, the extent of the problem is 500 km, the maximum height of the loads is 
   10000 m, and the number of load intervals is 120. 

4. Define heights and densities of loads. You can input these in a window or alternatively
   load them from a txt file (tab delimited) with two columns: 

   column 1: heights
   column 2: densities

   sample.txt is one of such files

5. Plot undeformed topography

6. Solve the problem. 

7. Plot deformed topography and deflection profile

8. Check the results: maximum deflection, distance to the forebulge, amplitude of forebulge,
   etc.

NOTE: When using an input file, keep track of the extent of the loaded region, the extent 
of the analysis, the maximum height of the loads, and the number of intervals used to 
discretize the load. You will need this info in step 3. The extent of the analysis can be 
equal or greater than the extent of the loaded region. The extent of the analysis cannot 
be lower than the extent of the loaded region.


RUNNING FLEX2D

In the command matlab window simply change to the flex2d folder

Type flex2d to run the program

Any problems email nfcd@mac.com

These scripts are copyright of Nestor Cardozo 2009. They can only be used for academic/research purposes and they cannot be distributed by third parties

Nestor Cardozo assumes no liability for damages, direct or consequential, which may result from the use of the scripts