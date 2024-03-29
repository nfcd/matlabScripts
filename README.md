# matlabScripts
A collection of Matlab scripts for fault related folding, basin analysis, and the book Structural Geology Algorithms: Vectors and Tensors (Allmendinger et al., 2012).

## Fault related folding

- faultBendFolding: 2D simple step, mode I fault bend folding
- constantThickness: 2D simple step, constant thickness fault propagation folding
- fixedAxis: 2D simple step, fixed axis fault propagation folding
- forward2D: 2D forward trishear fault propagation folding
- forwardP3D: Pseudo 3D forward trishear fault propagation folding
- forwardT3D: True 3D forward trishear fault propagation folding
- trishearCircFaultIS: 2D listric thrust. Trishear in front of thrust tip and inclined shear in backlimb
- trishearCircFaultPar: 2D listric thrust. Trishear in front of thrust tip and parallel shear in backlimb
- 3DtrishearCircFaultIS: 3D listric thrust. Trishear in front of thrust tip and inclined shear in backlimb.
- inverseBed: 2D, optimized trishear inverse modeling of folded beds (curved lines)
- inverseProf: 2D, optimized trishear inverse modeling of bedding intersections along a topographic profile
- inverseP3D: Pseudo 3D, optimized trishear inverse modeling of folded beds (folded surfaces)
- inverseT3D: True 3D, optimized trishear inverse modeling of folded beds (folded surfaces)

NOTE: You need the [Matlab Optimization Toolbox](http://www.mathworks.com/products/optimization/) to run the optimized trishear inversion scripts above.

The optimized trishear inversions above are highly affected by local minima (Cardozo and Aanonsen, 2009). To avoid being trapped in local minima, genetic and direct search algorithms can be used. The scripts below are new versions of 2D and 3D optimized trishear inversions using pattern search and simulated annealing algorithms (Cardozo et al. 2011).

- inverseBedPS: 2D, optimized trishear inverse modeling of folded bed (curved lines) using pattern search
- inverseBedSA: 2D, optimized trishear inverse modeling of folded bed (curved lines) using simulated annealing
- inverseProfPS: 2D, optimized trishear inverse modeling of bedding intersections along a topographic profile using pattern search
- inverseProfSA: 2D, optimized trishear inverse modeling of bedding intersections along a topographic profile using simulated annealing
- inverseP3DPS: Pseudo 3D, optimized trishear inverse modeling of folded beds (folded surfaces) using pattern search
- inverseP3DSA: Pseudo 3D, optimized trishear inverse modeling of folded beds (folded surfaces) using simulated annealing
- inverseT3DPS: True 3D, optimized trishear inverse modeling of folded beds (folded surfaces) using pattern search
- inverseT3DSA: True 3D, optimized trishear inverse modeling of folded beds (folded surfaces) using simulated annealing

NOTE: You need the [Matlab Global Optimization toolbox](https://se.mathworks.com/products/global-optimization.html) to run the optimized trishear inversions with pattern search or simulated annealing.

## Basin analysis:

- backstrip: 1D Airy backstripping
- flex2d: 2D flexural modeling. Continuous and broken beam of constant elastic thickness
- flex3d: 3D flexural modeling. Continuous plate of constant elastic thickness
- flex3dv: 3D flexural modeling. Continuous plate of variable elastic thickness (finite difference solution)

## Structural Geology Algorithms:

- StructuralGeologyAlgorithms: Matlab functions for the book Structural Geology Algorithms: Vectors and Tensors.

Any questions, please contact me [nestor.cardozo@uis.no](mailto:nestor.cardozo@uis.no).

