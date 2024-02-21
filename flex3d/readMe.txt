Flex3D solves for the deflection of an infinite, constant elastic thickness lithosphere, overlying a fluid astenosphere. The analogy is an infinite plate on an elastic foundation. The solution algorithm is that of Timoshenko and Woinowsky-Krieger (1959).

The workflow of the program is:

1. Start the program by typing flex3d in the Matlab command window. This will bring a figure window with title FLEX3D. It is in this window that you will interact with the program

2. Bring an image using the Image -> Load Image menu. The image is usually a map with contours of tectonic and sedimentary loads. Image testim.jpg is provided as an example.

3. Scale the image using the Image -> Scale Image menu. You will need to define the upper left and lower right corner of the image, and its horizontal and vertical dimensions.

4. Set up the plate geometry using the Parameters -> Set Plate Geometry menu. The plate is divided into a central region and four outer regions. The central region is the area of interest. It will cover the imported image. Outer regions are designed to avoid side effects. The user must define the number of columns and rows for the central and outer regions. For the central region, the density of the mesh should be based on the data and the resolution at which loads should be entered. For the outer regions, a bigger discretization interval can be used. 

5. Define elastic properties of the lithosphere using the Parameters -> Elastic Properties menu.

6. Define gravity and density of the astenosphere using the Parameters -> Gravity and Restoring Force menu.

6. Plot the plate and the image using the Plot -> Plate and Image menu. Zoom in to the area of interest (using the Matlab zoom tools) where you would like to input loads.

7. Interactively input loads using the Loads -> Load Input menu. Notice that for each cell in the plate where there are loads you need to click the cell and enter the height and density of the load. This can be cumbersome, and it can be better to modify the program such that loads can be input from a text file. I let that to the user. Notice that if you make mistakes entering loads, the loads can be removed using the Loads -> Load Removal menu. After that, just click the loads that need to be removed. Press the escape key to exit either the Input or Removal Load menus.

8. Plot the undeformed topography using the Plot -> Undeformed topography menu

9. Solve using the Solve -> Classical theory menu. Warning: Solving takes some time depending on the number of loads and cells in the plate. The analytical solution used in the program is slow. BE PATIENT. After solving the program will write in the command window the time it used to solve the problem. This will be the signal to continue with the next steps.

10. Plot the deformed topography using the Plot -> Deformed Topography menu. This will plot the deformed topography and an additional figure with the displacement in 3D.

11. Plot contours of deformed topography using the Plot -> Topographic contours menu. You will need to enter the desired contour values, as well as the x and y limits of the plot. The last ones should be obviously coincident with your area of interest. After that you can click the contours to label them, or just press enter.

12. Save loads using the Data -> Save Data menu. This will make a matlab file with the image, plate, and loads for posterior runs. Saved loads can be opened using the Data -> Load Data menu. After that you can just modify the loads using the Loads menus, and solve for the deflection.

The same workflow can be followed without an image in the background. Start the program, make your plate, define elastic properties and body forces, load it, solve for the loads, and plot.

Any questions, please email me at nfcd@mac.com

These scripts are copyright of Nestor Cardozo 2009. They can only be used for academic/research purposes and they cannot be distributed by third parties

Nestor Cardozo assumes no liability for damages, direct or consequential, which may result from the use of the scripts

