
# Tutorial: End-to-end analysis of the Shallow Water Equations with LSPG: Writing the LSPG ROM 
This page walks through writing a coupler that interfaeces with the LSPG.

The first step in our LSPG driver file is to define and initialize the FOM. To do this, we define the domain, number of grid points, load in parameter, etc. 
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_rom/run_lspg.cc
68:90
```
After extracting the relevant types, we now define a parameter grid on which we will solve the FOM, and then initialize the app. The interface required to initialize the app can be viewed in the swe2d.hpp source file, and requires 
arguments for the lengnth of the domain and number of grid points in the x and y directions, as well as an array that sets the system parameters. 
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/offline_phase/run_fom_for_training_params.cc
37:51
```
Next, we construct a Crank Nicolson time stepper that we will use to march the problem in time. In Pressio, the steppers (1) act on Pressio datatypes that wrap the native datatype and (2) are templated on the time scheme types for the state, residual, jacobian, and application. As such, we first define the relevant Pressio wrapped data types, and then define the stepper type. 
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/offline_phase/run_fom_for_training_params.cc
54:62
```
Next, we define the linear solver type, and construct the solver. Here, we use the stabilized biconjugate gradient method. 
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/offline_phase/run_fom_for_training_params.cc
64:67
```
Finally, we define the relevant information for our time grid, loop over the parameter instances, and then solve the FOM.
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/offline_phase/run_fom_for_training_params.cc
68:103
```


