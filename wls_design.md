# Required Functions for WLS

variables:

 - nts: width of time stencil, for forward euler nts = 2
 - ns:  number of time steps in the window
 - K: ROM basis size
 - N: FOM size

 
### Evauate single residual
This computes the residual at a single step, and is the same as LSPG.
Required inputs are:

 1. The nts states required for the time discretization, type = rom vectors
 2. The time instance (for non-autonomous functions), type=double
 3. The time step, type=double

Outputs/required storage is

 1. residual, FOM vector of size N  

### Evauate windowed residual
This computes the residual over the time window. Presumably this will use the "evaluate single residual" function. 
Required inputs are:

 1. The number of time steps, ns, to compute the residual over, type=int
 2. The nts + ns - 1 states required for the time discretization, type = rom vectors
 3. The time instances (for non-autonomous functions), type=vector of doubles
 4. The time steps, type = vector of doubles

Outputs/required storage is

 1. windowed residual, extended FOM vector of size N * ns

### Evauate single Jacobian
This computes the Jacobian of a single step, and is the same as LSPG
Required inputs are:

 1. The nts states required for the time discretization, type = rom vectors
 2. The time instance (for non-autonomous functions), type=double
 3. The time step, type=double
 4. The argument of the residual with which to compute the Jacobian of, type=int

Outputs/required storage is

 1. Jacobian, FOM matrix of size N * K

### Evauate windowed Jacobian
This computes the Jacobian of over the time window. Presumably this will use the "evaluate single Jacobian" function. 
Required inputs are:

 1. The number of time steps, ns, to compute the residual over, type=int
 2. The nts + ns - 1 states required for the time discretization, type = rom vectors
 2. The time instance (for non-autonomous functions), type=double
 3. The time step, type=double
 4. The argument of the residual with which to compute the Jacobian of, type=int

 Outputs/required storage is

 1. Jacobian, extended FOM sparse matrix of size N ns* K ns