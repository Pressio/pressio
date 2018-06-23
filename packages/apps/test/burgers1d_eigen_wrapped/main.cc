
#include "ode_observer.hpp"

#include "vector/core_vector_serial_eigen.hpp"
#include "matrix/core_matrix_dense_serial_eigen.hpp"
// app class
#include "apps_burgers1d_eigen.hpp"
// integrator
#include "integrators/ode_integrate_n_steps.hpp"
// steppers
#include "step_methods/ode_explicit_euler_stepper.hpp"
#include "step_methods/ode_explicit_runge_kutta4_stepper.hpp"
#include "step_methods/ode_implicit_euler_stepper.hpp"
// solvers
#include "experimental/solvers_linear_eigen.hpp"
#include "experimental/solvers_nonlinear_newton_raphson.hpp"


struct stateResizer{
  using native_state_t = apps::burgers1dEigen::state_type;
  using vec_t = core::vector<native_state_t>;
  void operator()(const vec_t & source, vec_t & dest){
    if ( dest.size()!=source.size() )
      dest.resize(source.size());
  };
};


int main(int argc, char *argv[])
{
  using native_state_t = apps::burgers1dEigen::state_type;
  using native_jac_t = apps::burgers1dEigen::jacobian_type;
  using scalar_t = apps::burgers1dEigen::scalar_type;
  using model_eval_t = apps::burgers1dEigen;

  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  model_eval_t appObj(mu, 10);
  appObj.setup();

  // integrate in time starting from y0
  scalar_t final_time = 35.;
  scalar_t dt = 0.01;

  using state_t = core::vector<native_state_t>;
  native_state_t y0n = appObj.getInitialState();
  // {
  //   //********************************************
  //   // EXPLICIT EULER
  //   //********************************************
  //   state_t y0(y0n);
  //   snapshot_collector collectorObj;

  //   using stepper_t = ode::explicitEulerStepper<
  //     state_t, state_t, scalar_t, stateResizer,
  //     model_eval_t, scalar_t /*default = standard policy */>;
  //   stepper_t stepperObj(appObj);
    
  //   ode::integrateNSteps(stepperObj, y0, collectorObj,
  // 			 0.0, dt, final_time/dt);
  //   std::cout << collectorObj.getCount() << std::endl;

  //   std::cout << "Final solution " << std::endl;
  //   for (int i=0; i<y0.size(); ++i)
  //     std::cout << std::setprecision(14) << y0[i]  << " ";
  //   std::cout << std::endl;
  // }
  // {
  //   //********************************************
  //   // EXPLICIT RUNGE KUTTA 4
  //   //********************************************
  //   state_t y0(y0n);
  //   snapshot_collector collectorObj;

  //   using stepper_t = ode::explicitRungeKutta4Stepper<
  //     state_t, state_t, scalar_t, stateResizer,
  //     model_eval_t, scalar_t /*default = standard policy */>;
  //   stepper_t stepperObj(appObj);

  //   ode::integrateNSteps( stepperObj, y0, collectorObj,
  // 			  0.0, dt, final_time/dt );
  //   std::cout << collectorObj.getCount() << std::endl;

  //   std::cout << "Final solution " << std::endl;
  //   for (int i=0; i<y0.size(); ++i)
  //     std::cout << std::setprecision(14) << y0[i]  << " ";
  //   std::cout << std::endl;
  // }

  //  {
    using jac_t = core::matrix<native_jac_t>;

    //********************************************
    // IMPLICIT EULER
    //********************************************
    state_t y0(y0n);
    snapshot_collector collectorObj;

    // linear solver
    using lin_solve_t =
      solvers::experimental::linearSolver<jac_t,state_t,state_t>;
    lin_solve_t ls;
    
    // nonlinear solver
    using nonlin_solve_t =
      solvers::experimental::newtonRaphson<state_t,state_t,jac_t,lin_solve_t>;
    nonlin_solve_t nonls(ls);

    //stepper
    using stepper_t = ode::implicitEulerStepper<
      state_t, state_t, jac_t, scalar_t, stateResizer,
      model_eval_t, scalar_t, nonlin_solve_t
      /*, default = standard policy for res and jac*/>;
    stepper_t stepperObj(appObj, nonls);

    ode::integrateNSteps( stepperObj, y0, 0.0, dt, final_time/dt );
    // std::cout << collectorObj.getCount() << std::endl;

    // // std::cout << "Final solution " << std::endl;
    for (int i=0; i<y0.size(); ++i)
      std::cout << std::setprecision(14) << y0[i]  << " ";
    std::cout << std::endl;
    // //  }

  
  return 0;
}
