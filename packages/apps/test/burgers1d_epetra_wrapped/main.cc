
#include "CORE_ALL"
#include "ODE_ALL"
//#include "SOLVERS_EXP"
// app class
#include "apps_burgers1d_epetra.hpp"
#include "ode_observer.hpp"
//#include "../apps_helper_ode.hpp"

#include "experimental/solvers_l2_vector_norm.hpp"
#include "experimental/solvers_linear_base.hpp"
#include "experimental/solvers_linear_dense.hpp"
#include "experimental/solvers_linear_factory.hpp"
#include "experimental/solvers_linear_iterative_policies_trilinos.hpp"
#include "experimental/solvers_linear_iterative_traits.hpp"
#include "experimental/solvers_linear_iterative.hpp"
#include "experimental/solvers_nonlinear_base.hpp"
#include "experimental/solvers_nonlinear_factory.hpp"
#include "experimental/solvers_nonlinear_iterative.hpp"
#include "experimental/solvers_nonlinear_traits.hpp"
#include "experimental/solvers_policy_linear_dense_eigen.hpp"
#include "experimental/solvers_policy_linear_iterative_eigen.hpp"
#include "experimental/solvers_policy_linear_iterative_trilinos.hpp"
#include "experimental/solvers_policy_nonlinear_iterative.hpp"



int main(int argc, char *argv[])
{
  //-------------------------------
  // define native types
  using native_state_t = apps::Burgers1dEpetra::state_type;
  using native_space_residual_t = apps::Burgers1dEpetra::space_residual_type;
  using native_jac_t = apps::Burgers1dEpetra::jacobian_type;
  using scalar_t = apps::Burgers1dEpetra::scalar_type;
  using target_app_t = apps::Burgers1dEpetra;
  
  // MPI_Init(&argc,&argv);
  // int rank; // My process ID
  // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // Epetra_MpiComm Comm(MPI_COMM_WORLD);

  // //-------------------------------
  // // create app object
  // std::vector<double> mu({5.0, 0.02, 0.02});
  // target_app_t appObj(mu, 20, &Comm);
  // appObj.setup();
  // auto & y0n = appObj.getInitialState();
  // auto & r0n = appObj.getInitialResidual();
  // auto & j0n = appObj.getInitialJacobian();

  {
    // //-------------------------------
    // // wrapper types 
    // using state_t = core::Vector<native_state_t>;
    // using jac_t = core::Matrix<native_jac_t>;
    // using res_t = core::Vector<native_space_residual_t>;

    // state_t y(y0n);
    // space_res_t r(r0n);
    // jac_t j0FOM(j0n);
    
    // // stepper
    // using stepper_t = ode::ImplicitEulerStepper<
    //   state_t, res_t, jac_t, target_app_t>;
    // // using stepper_t = ode::ExplicitRungeKutta4Stepper<
    // //   state_t, space_res_t, scalar_t, target_app_t>;
    // stepper_t stepperObj(appObj, y, r);

    // Create linear solver using a valid solver type
    auto solver = solvers::NonLinearSolvers::createSolver<solvers::nonlinear::NewtonRaphson>();
    
    // // integrate in time 
    // snapshot_collector collObj;
    // scalar_t final_t = 35;
    // scalar_t dt = 0.01;
    // ode::integrateNSteps(stepperObj, y, 0.0, dt, static_cast<unsigned int>(dt/dt), solver);
    // //y.data()->Print(std::cout << std::setprecision(14));
    // // // //printSol("Exp Euler", y0);
  }


  // {
  //   //-------------------------------
  //   // wrapper types 
  //   using state_t = core::Vector<native_state_t>;
  //   using jac_t = core::Matrix<native_jac_t>;
  //   using space_res_t = core::Vector<native_space_residual_t>;

  //   state_t y(y0n);
  //   space_res_t r(r0n);

  //   // stepper
  //   using stepper_t = ode::ExplicitEulerStepper<
  //     state_t, space_res_t, target_app_t>;
  //   // using stepper_t = ode::ExplicitRungeKutta4Stepper<
  //   //   state_t, space_res_t, scalar_t, target_app_t>;
  //   stepper_t stepperObj(appObj, y, r);

  //   // integrate in time 
  //   snapshot_collector collObj;
  //   scalar_t final_t = 35;
  //   scalar_t dt = 0.01;
  //   ode::integrateNSteps(stepperObj, y, 0.0, dt, static_cast<unsigned int>(final_t/dt), collObj);
  //   y.data()->Print(std::cout << std::setprecision(14));
  //   // // //printSol("Exp Euler", y0);
  // }
  
  MPI_Finalize();
  return 0;
}
