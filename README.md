# Designing Control Barrier Functions Using a Dynamic Backup Policy
This repository contains MATLAB source code for the paper [Designing Control Barrier Functions
Using a Dynamic Backup Policy].


For more information about our work, please visit [ROCC Team@CU Boulder](https://www.colorado.edu/faculty/nicotra/robotics-optimization-and-constrained-control).

## MATLAB
### Dependencies & Installation
The following MATLAB toolboxes are required:
* [Optimization Toolbox](https://www.mathworks.com/products/optimization.html)
* [Control Systems Toolbox](https://www.mathworks.com/products/control.html)

The following tools are recommended:
* [MOSEK](https://www.mosek.com/) as a specialized convex QP solver. Free [academic licenses](https://www.mosek.com/products/academic-licenses/) are available.

The package is lightweight and there is no installation beyond adding the folder to
your path:
```
addpath(userpath+"\traj-DSM-CBFs")
```

To run the code, execute the file `pend_sim_all.m`.

## Acknowledgements
This work was supported by the University of Colorado Boulder and the NSF-CMMI Award #2411667.
