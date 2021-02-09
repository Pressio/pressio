
# Tutorial: End-to-end analysis of the Shallow Water Equations with LSPG 

@m_class{m-block m-info}

@par
This tutorial works through an end-to-end analysis where we use the LSPG approach to accelerate a forward model of the shallow water equations (SWEs). In this tutorial, we will:
1. Interface a SWE solver, written with Eigen data structures, to Pressio.
2. Use Pressio's time marching schemes to execute solves of the SWEs to construct training data.
3. Use supporting Python scripts to analyze the training data and construct ROM basis vectors.
4. Construct and run a standard LSPG ROM for novel parameter training instances.
5. Construct and run a hyperreduced LSPG ROM for novel parameter training instances.

@image html swetut_f1.gif width=90%

# The shallow water equations
We consider the shallow water equations on the spatial domain @f$\Omega = [-\frac{L}{2},\frac{L}{2}] \times  [-\frac{L}{2},\frac{L}{2}]@f$, which comprises the PDE system
@f[
\begin{split}
&\frac{\partial h}{\partial t} + \frac{\partial}{\partial x }(  h u) + \frac{\partial}{\partial y }( h v) = 0\\
&\frac{\partial h u}{\partial t} + \frac{\partial}{\partial x} (h u^2 + \frac{1}{2} \mu_1 h^2) + \frac{\partial}{\partial y }( h u v) = \mu_3 hv\\
&\frac{\partial h v}{\partial t} + \frac{\partial}{\partial x} (h u v) + \frac{\partial}{\partial y }( h v^2 +  \frac{1}{2} \mu_1 h^2) = \mu_3 hu.
\end{split}
@f]
In the above, @f$h : \Omega \rightarrow \mathbb{R}@f$ is the height of the water surface,  @f$u : \Omega \rightarrow \mathbb{R}@f$ is the x-velocity, and @f$v : \Omega \rightarrow \mathbb{R}@f$ is the y-velocity. The system has three parameters: 
@f$\mu_1 \in [3,9]@f$ is the gravity parameter, @f$\mu_2@f$ controls the magnitude of the initial pulse, and @f$\mu_3@f$ controls the magnitude of the Coriolis forcing. 

# Offline phase
In the offline phase, we
1. Interface with and run the full-order model
2. Extract training data and perform proper orthogonal decomposition to find a basis for the ROM
3. Select indices for hyper-reduction by Q-sampling
## The full-order model 
A full-order model that can be used to solve the SWEs is located in the packages/apps section of the Pressio repo (see [here](https://github.com/Pressio/pressio/tree/swe2d/packages/apps/src/swe2d) ). The full-order model employs a first-order finite volume 
discretization with the Rusanov flux scheme at the cell interfaces.  

The first step in our analysis is to run the full-order model for training parameter instances. To do this, we can write the driver file *run_fom_for_training_params.cc*. Thie file couples Pressio to the application, and uses Pressio's time marching scheme to solve the model. See [here](./md_pages_tutorials_tutorial3_fom.html) for a step-by-step walkthrough of constructing this driver file.

To run the driver file, move to the offline_phase directory and run the script: 
```bash
cd tutorials/swe2d/offline_phase
./run_fom_for_training_params
```

To do this, we are going to write a driver file that couples Pressio to the application, and use Pressio's time-marching schemes to solver the model. 

The first step in our driver file is to define types. Here, we start by defining our application to be the swe2d app, and then extract the relavant types from the app.
```cpp
@codesnippet
../../../../pressio-tutorials/tutorials/swe2d/offline_phase/run_fom_for_training_params.cc
30:37
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

# How to create a default Galerkin problem?
To create a default Galerkin problem object, one needs:
1. a FOM object satisfying the API described [here]()
2. a linear decoder (see [this tutorial](./md_pages_tutorials_tutorial1.html))
3. a rom state
4. a FOM reference state
5. [optional] an object with specific kernels when the FOM types are not natively supported by pressio.<br/>
See [this page](./md_pages_custom_ops_default_gal_exp.html) for more information on this.

Synopsis:

```cpp
using ode_tag = pressio::ode::explicitmethods::Euler;

using pressio::rom::galerkin::createDefaultProblem;
auto Problem = createDefaultProblem<ode_tag>(fomObj, decoder, yRom, yRef, [, opsObject]);
```
Note the function is templated on the ode tag.
To select a different time stepping scheme, one can change the tag.
To see the list of currently supported explicit stepping schemes, see \todo.


# How to solve a default Galerkin problem?

Once the target problem object is created, the reduced system
can be integrated in time. Here we provide the most basic function
to do so, which advances the system for a fixed number of steps.
Synopsis:

```cpp
// solve for fixed number of steps and time step
pressio::rom::galerkin::solveNSteps(problem,     # problem object
								    yRom,        # rom state to advance
								    t0,          # initial time
									dt,          # time step
									Nsteps       # number of steps
									[, observer] # optional observer
								   )
```
The optional argument allows one to pass an "observer" object whose
purpose is to monitor the evolution of the reduced state.
The observer is called back by pressio4py during the time integration
at every time step. This can be useful to, e.g., save the
generalized coordinates, or usign them to perfom some other operation.


# Putting all steps together

```cpp
// create adapter/FOM object
// ...

// create the decoder
decoder = /*see, e.g., tutorial for linear decoder */

// create ROM state: here we use Eigen
using rom_state_t = pressio::containers::Vector<Eigen::VectorXd>;
rom_state_t yRom(/*whatever rom size needed*/);

// create the Galerkin problem
using ode_tag = pressio::ode::explicitmethods::Euler;

using pressio::rom::galerkin::createDefaultProblem;
auto Problem = createDefaultProblem<ode_tag>(fomObj, decoder, yRom, yRef, [, opsObject]);

// to solve, we set as an example, t0=0, dt=0.1, Nstep = 100
pressio::rom::galerkin::solveNSteps(Problem, yRom, 0., 0.1, 100)
```
Note the function is templated on the ode tag.
To select a different time stepping scheme, one can change the tag.
To see the list of currently supported explicit stepping schemes, see \todo.



<!-- The observer class must meee the following API: -->
<!-- ```py -->
<!-- class OdeObserver: -->
<!--   def __init__(self): pass -->

<!--   def __call__(self, timeStep, time, romState): -->
<!-- 	# do what you want with romState -->
<!-- ``` -->
<!-- Note that we are working on enriching the API to integrate in time. -->
<!-- For example, we will soon support function class to advance the problem -->
<!-- until a condition is met, or until a target time is reached. -->


<!-- # Want to see all the above pieces in action? -->

<!-- Look at [this demo](./md_pages_demos_demo1.html) that uses -->
<!-- default Galerkin for a 1d PDE. -->


<!-- # Some considerations -->
<!-- @m_class{m-block m-warning} -->

<!-- @par -->
<!-- One might wonder how the above formulation can be efficient, -->
<!-- given that the right-hand side of the reduced system scales -->
<!-- with the FOM degrees of freedom. -->
<!-- This is true: the reduced system obtained from a -->
<!-- *default* problem reduces the spatial degrees of freedom, -->
<!-- but is typically not efficient because at every evaluation of the RHS, -->
<!-- it requires a large matrix vector product. -->
<!-- Thus, a default Galerkin is typically used for exploratory -->
<!-- analysis when computational efficiency is **not** a primary -->
<!-- goal, e.g. to test the feasibility of ROMs for a target problem, -->
<!-- or try different basis. -->
<!-- When computational efficiency is critical, one needs to -->
<!-- resort to hyper-reduction techniques to reduce the cost of the matrix-vector -->
<!-- product. This is covered in subsequent tutorials. -->
