
# Tutorial: End-to-end analysis of the Shallow Water Equations with LSPG 
The first step in our driver file is to define types. Here, we start by defining our application to be the swe2d app, and then extract the relavant types from the app.
```cpp
  // Define types to use
  using app_t	        = ::pressio::apps::swe2d;
  using scalar_t        = double;
  using app_state_t	= typename app_t::state_type;
  using app_rhs_t	= typename app_t::velocity_type;
  using app_jacob_t	= typename app_t::jacobian_type;

  // Create arrays for parameter sweep
```
After extracting the relevant types, we now define a parameter grid on which we will solve the FOM, and then initialize the app. The interface required to initialize the app can be viewed in the swe2d.hpp source file, and requires 
arguments for the lengnth of the domain and number of grid points in the x and y directions, as well as an array that sets the system parameters. 
```cpp
  // Create arrays for parameter sweep
  std::vector<double> params0Array{3.,6.,9.}; //parameter controls gravity
  double              params1 = 0.125; //parameter controls the magnitude of the initial pulse (kept constant here)
  std::vector<double> params2Array{0.05,0.15,0.25}; // parameter controls the forcing

  // Construct the app
  int nx = 128;
  int ny = 128;
  scalar_t Lx = 5;
  scalar_t Ly = 5;
  scalar_t params[3];
  params[0] = params0Array[0]; 
  params[1] = params1; 
  params[2] = params2Array[0]; 
  app_t appObj(Lx,Ly,nx,ny,params);
```
Next, we construct a Crank Nicolson time stepper that we will use to march the problem in time. In Pressio, the steppers (1) act on Pressio datatypes that wrap the native datatype and (2) are templated on the time scheme types for the state, residual, jacobian, and application. As such, we first define the relevant Pressio wrapped data types, and then define the stepper type. 
```cpp
  using ode_state_t = pressio::containers::Vector<app_state_t>;
  using ode_res_t   = pressio::containers::Vector<app_rhs_t>;
  using ode_jac_t   = pressio::containers::SparseMatrix<app_jacob_t>;
  using ode_tag = pressio::ode::implicitmethods::CrankNicolson;

  // define type for stepper object
  using stepper_t = pressio::ode::ImplicitStepper<
    ode_tag, ode_state_t, ode_res_t, ode_jac_t, app_t>;

```
Next, we define the linear solver type, and construct the solver. Here, we use the stabilized biconjugate gradient method. 
```cpp
  using lin_solver_t = pressio::solvers::linear::Solver<
    pressio::solvers::linear::iterative::Bicgstab, ode_jac_t>;
  lin_solver_t linSolverObj;

```
Finally, we define the relevant information for our time grid, loop over the parameter instances, and then solve the FOM.
```cpp
  // Define information for time stepping
  scalar_t t = 0; 
  scalar_t et = 10.; 
  scalar_t dt = 0.02;
  auto Nsteps = static_cast<::pressio::ode::types::step_t>(et/dt);
  int fileNo = 0;
  for (int params0Counter = 0; params0Counter < params0Array.size(); params0Counter++){
    for (int params2Counter = 0; params2Counter < params2Array.size(); params2Counter++){
      // Set parameters
      params[0] = params0Array[params0Counter];
      params[2] = params2Array[params2Counter];
      appObj.setParams(params);

      // Construct initial state
      ode_state_t y(appObj.getGaussianIC(params[1]));

      // Construct stepper
      stepper_t stepperObj(y, appObj);

      // Construct observer that saves output
      std::string filename = "solution";
      filename += std::to_string(fileNo);
      filename += ".bin";
      observer<ode_state_t> Obs(filename);

      // Create nonlinear solver
      auto NonLinSolver=
        pressio::solvers::nonlinear::createNewtonRaphson(stepperObj, y,linSolverObj);
      NonLinSolver.setTolerance(1e-11);

      // Solve problem
      pressio::ode::advanceNSteps(stepperObj, y, t, dt, Nsteps,Obs, NonLinSolver);
      Obs.closeFile();
      fileNo += 1;
    }
  }
```


