
# Tutorial: End-to-end analysis of the Shallow Water Equations with LSPG: Writing the LSPG ROM 
This page walks through writing a driver file for an LSPG ROM of the shallow water equations. The full code for this coupler is available in the [pressio-tutorials repo]((https://github.com/Pressio/pressio-tutorials/blob/swe2d_tutorial/tutorials/swe2d/online_phase/lspg_rom/run_lspg.cc). The first step in writing our driver is to include the relevant headers from the Pressio library. Here, we include [pressio_apps.hpp](https://github.com/Pressio/pressio/blob/master/packages/pressio_apps.hpp), and [pressio_rom_lspg.hpp](https://github.com/Pressio/pressio/blob/master/packages/pressio_lspg.hpp), which enable us to access our Shallow Water application and LSPG capabilities, respectively. 
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_rom/run_lspg.cc
1:2
```

The next step for us is to start out int main file, and first we define and initialize the FOM. To do this, we define the domain, number of grid points, load in parameter, etc. 

```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_rom/run_lspg.cc
68:87
```

After constructing the app object, we now read in the bases. For standard LSPG, we need to read in the bases for all mesh points. Here, we assume the basis exists in a local file named *basis.txt*. The following lines of code read in the basis:

```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_rom/run_lspg.cc
89:97
```
We note that the above code uses the *convertFromVVecToMultiVec* and *readBasis* functions, which are defined at the top of the script.
 
Now that we have read in the basis, we create a Pressio decoder object. In this instance where we have a linear basis, the decoder, in essence, computes the product @f$\boldsymbol \Phi \hat{\boldsymbol x} @f$. We construct a decoder as follows
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_rom/run_lspg.cc
98:104
```
Next, we set the reference state for the ROM, @f$\boldsymbol x^{ref}@f$. Note that this reference state will be used when we reconstruct the state vector, i.e., the state will be reconstructed as
@f[\hat{ \boldsymbol x} = \boldsymbol \Phi \hat{\boldsymbol x} + \boldsymbol x^{ref}@f]
In this case we set the reference state to be the initial condition, and set the reference states as
 ```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_rom/run_lspg.cc
106:109
```
We now begin constructing the ROM problem. We start by setting our LSPG state type to be an eigen vector. We note that this is decoupled from the FOM, and we could alternatively use, e.g., Kokkos. In the following snippet of code, we set this type, initialize the ROM state, and set it to be zero. Note that this last step is due to the fact that we employ the initial condition as our ROM reference state, so the initial conditions for our ROM state are zero.
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_rom/run_lspg.cc
111:118
```

We now construct the LSPG problem. To do this, we have to define the time marching scheme, which here we set to be the second-order Crank Nicolson method. We define this time stepping scheme and construct the LSPG ROM as follows
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_rom/run_lspg.cc
120:123
```

Now, we need to define the linear and nonlinear solvers use in our LSPG problem. Here, we employ a least-squares conjugate gradient solver, which is constructed via
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_rom/run_lspg.cc
125:130
```
With these elements defined, we can define a Gauss--Newton solver for our LSPG ROM. We do this in the following snippet off code, and additionally set the convergence tolerance and maximum number of iterations
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_rom/run_lspg.cc
132:135
```
Next we define information for our time marching. We set the initial time, time step, end time, and total number of time steps, respectively,
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_rom/run_lspg.cc
137:140
```

Next, we construct an observer. Observers act like hooks, and are called at the end of each time step. The observer itself is defined outside of int main, and here looks like:
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_rom/run_lspg.cc
34:61
```
When initialized, our observer writes the reference state to file in *state_ref.bin*. When this observer is called after each time step, it will then write the ROM solution to a *solution.bin* file. We initialize this object inside our int main as follows:

```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_rom/run_lspg.cc
142:143
```

Finally, we can now use Pressio to solve the ROM. The following code advances the ROM in time for *Nsteps* time steps. Note that, as LSPG consists of solving a series of residual minimization problems, we call this *solveNSequentialMinimizations*:
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/online_phase/lspg_rom/run_lspg.cc
145:150
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

This completes our description of writing the LSPG coupler. Click [here](./md_pages_tutorials_tutorial3.html) to return to the SWE tutorial.


