
# Tutorial: End-to-end analysis of the Shallow Water Equations with LSPG: Coupling to the FOM 
This page walks through writing a driver file for a FOM of the shallow water equations using Pressio's time marching schemes. The full code for this coupler is available in the [pressio-tutorials repo](https://github.com/Pressio/pressio-tutorials/blob/swe2d_tutorial/tutorials/swe2d/offline_phase/run_fom_for_training_params.cc).
The first step in our driver file is to define types. Here, we start by defining our application to be the swe2d app, and then extract the relevant types from the app.
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/offline_phase/run_fom_for_training_params.cc
30:37
```
After extracting the relevant types, we now define a parameter grid on which we will solve the FOM, and then initialize the app. The interface required to initialize the app can be viewed in the swe2d.hpp source file, and requires 
arguments for the length of the domain and number of grid points in the x and y directions, as well as an array that sets the system parameters. 
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/offline_phase/run_fom_for_training_params.cc
37:51
```
Next, we construct a Crank Nicolson time stepper that we will use to march the problem in time. In Pressio, the steppers (1) act on Pressio data types that wrap the native datatype and (2) are templated on the time scheme types for the state, residual, Jacobian, and application. As such, we first define the relevant Pressio wrapped data types, and then define the stepper type. 
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/offline_phase/run_fom_for_training_params.cc
54:62
```
Next, we define the linear solver type, and construct the solver. Here, we use the stabilized bi-conjugate gradient method. 
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

This completes our description of writing the FOM coupler. Click [here](./md_pages_tutorials_tutorial3.html) to return to the SWE tutorial.


