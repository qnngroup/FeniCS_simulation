# FeniCS_simulation

Author: Adina Bechhofer (2021 - 2022)

### General Description 

This project implements a FEniCS simulation for a PVNCT emitter. We are also trying to get an adjoint method going for an efficient optimization of emission current. The idea is to treat the coefficients of the finite element basis as the simulation varables instaed of the voltages at the nodes themselves.  
This project depends on Python3, fenics, dolfin, numpy, matplotlib, and scipy.

### Running 

`spline_geom.py` constructs an emitter as a spline interpolation of 88 boundry points and runs the poisson solver with specified boundary conditions. Then, the scrip bcalculates the electric field everywhere and uses it to calculate the Fowler Nordheim emission.  
`adjoint_with_mesh_pert.py` constructs the dg/dp term in the adjoint formulation by perturbing the mesh and recalculating the finite element matrix.  

