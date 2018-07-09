
#include "ode_observer.hpp"

#include "vector/concrete/core_vector_serial_eigen.hpp"
#include "matrix/concrete/core_matrix_dense_serial_eigen.hpp"
#include "matrix/concrete/core_matrix_sparse_serial_eigen.hpp"
// #include "CORE_VECTOR"
// #include "CORE_MATRIX"
// app class
#include "apps_burgers1d_eigen.hpp"
// integrator
#include "integrators/ode_integrate_n_steps.hpp"
// policies
#include "policies/custom/ode_residual_increment_policy.hpp"
#include "policies/custom/ode_jacobian_increment_policy.hpp"
// steppers
#include "steppers/explicit_steppers/ode_explicit_euler_stepper.hpp"
#include "steppers/explicit_steppers/ode_explicit_runge_kutta4_stepper.hpp"
#include "steppers/implicit_steppers/ode_implicit_euler_stepper.hpp"
//solvers
#include "experimental/solvers_linear_eigen.hpp"
#include "experimental/solvers_nonlinear_newton_raphson.hpp"


struct mysizer{
 using state_t = core::vector<apps::burgers1dEigen::state_type>;
 static size_t getSize(const state_t & obj){
   return obj.size();
 };
 static void resize(state_t & obj, size_t newSize){
   obj.resize(newSize);
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

  // integrate in time startxbi5ng from y0
  scalar_t dt = 0.001;
  scalar_t final_t = 35;

  using state_t = core::vector<native_state_t>;
  using jac_t = core::matrix<native_jac_t>;
  using residual_t = state_t;
  native_state_t y0n = appObj.getInitialState();
  snapshot_collector collectorObj;

  state_t y0(y0n);

  using res_pol_t = ode::policy::incrementBasedResidual<
    state_t, residual_t, model_eval_t, scalar_t, mysizer>;
  res_pol_t resObj(y0);

  using jac_pol_t = ode::policy::incrementBasedJacobian<
    state_t, jac_t, model_eval_t, scalar_t, mysizer>;
  jac_pol_t jacObj(y0);    
  
  // {
  //   //********************************************
  //   // EXPLICIT EULER
  //   //********************************************
  //   using stepper_t = ode::explicitEulerStepper<
  //     state_t, residual_t, scalar_t,
  //     model_eval_t, scalar_t, mysizer, res_pol_t>;
  //   stepper_t stepperObj(appObj, resObj);
    
  //   ode::integrateNSteps(stepperObj, y0, 0.0, dt, final_t/dt, collectorObj);
  //   //    std::cout << collectorObj.getCount() << std::endl;
  //   std::cout << "Final solution " << std::endl;
  //   for (int i=0; i<y0.size(); ++i)
  //     std::cout << std::setprecision(14) << y0[i]  << " ";
  //   std::cout << std::endl;
  // }

  {
    *y0.data() = y0n;
    //********************************************
    // EXPLICIT RK4
    //********************************************
    using stepper_t = ode::explicitRungeKutta4Stepper<
      state_t, residual_t, scalar_t,
      model_eval_t, scalar_t, mysizer, res_pol_t>;
    stepper_t stepperObj(appObj, resObj);
    
    ode::integrateNSteps(stepperObj, y0, 0.0, dt, final_t/dt, collectorObj);
    //    std::cout << collectorObj.getCount() << std::endl;
    std::cout << "Final solution " << std::endl;
    for (int i=0; i<y0.size(); ++i)
      std::cout << std::setprecision(14) << y0[i]  << " ";
    std::cout << std::endl;
  }
  
  // {
  //   *y0.data() = y0n;
  //   //********************************************
  //   // EXPLICIT RK4
  //   //********************************************
  //   using stepper_t = ode::explicitRungeKutta4Stepper<
  //     state_t, residual_t, scalar_t,
  //     model_eval_t, scalar_t, mysizer>;
  //   stepper_t stepperObj(appObj);
    
  //   ode::integrateNSteps(stepperObj, y0, 0.0, dt, final_t/dt, collectorObj);
  //   //    std::cout << collectorObj.getCount() << std::endl;
  //   std::cout << "Final solution " << std::endl;
  //   for (int i=0; i<y0.size(); ++i)
  //     std::cout << std::setprecision(14) << y0[i]  << " ";
  //   std::cout << std::endl;
  // }

  // {
  //   *y0.data() = y0n;

  //   //********************************************
  //   // IMPLICIT EULER
  //   //********************************************
  //   // linear solver
  //   using lin_solve_t = solvers::experimental::linearSolver<jac_t,state_t,state_t>;
  //   lin_solve_t ls;
    
  //   // nonlinear solver
  //   using nonlin_solve_t =
  //     solvers::experimental::newtonRaphson<state_t,state_t,jac_t,lin_solve_t>;
  //   nonlin_solve_t nonls(ls);

  //   //stepper 1
  //   using stepper_t = ode::implicitEulerStepper<
  //     state_t, residual_t, jac_t, scalar_t,
  //     model_eval_t, scalar_t, mysizer, nonlin_solve_t, res_pol_t, jac_pol_t>;
  //   stepper_t stepperObj(appObj, nonls, resObj, jacObj);
    
  //   ode::integrateNSteps( stepperObj, y0, 0.0, dt, final_t/dt, collectorObj);
  //   // // // // std::cout << collectorObj.getCount() << std::endl;
  //   // std::cout << "Final solution " << std::endl;
  //   std::cout << "\n-----------" << std::endl;
  //   std::cout << *y0.data();
  //   // // // for (int i=0; i<y0.size(); ++i)
  //   // // //   std::cout << std::setprecision(14) << y0[i]  << " ";
  //   // // // std::cout << std::endl;
  // }

  // {
  //   using jac_t = core::matrix<native_jac_t>;

  //   // linear solver
  //   using lin_solve_t =
  //     solvers::experimental::linearSolver<jac_t,state_t,state_t>;
  //   lin_solve_t ls;
    
  //   // nonlinear solver
  //   using nonlin_solve_t =
  //     solvers::experimental::newtonRaphson<state_t,state_t,jac_t,lin_solve_t>;
  //   nonlin_solve_t nonls(ls);

  //   // first define the auxiliary stepper for initial stepping
  //   // using aux_res_pol_t = ode::policy::implicitEulerIncrementResidual<
  //   //   state_t, residual_t, model_eval_t, scalar_t>;
  //   // aux_res_pol_t auxResObj(y0);
  //   // using aux_jac_pol_t = ode::policy::implicitEulerIncrementJacobian<
  //   //   state_t, jac_t, model_eval_t, scalar_t>;
  //   // aux_jac_pol_t auxJacObj(y0);
  //   // //stepper
  //   // using aux_stepper_t = ode::implicitEulerStepper<
  //   //   state_t, residual_t, jac_t, scalar_t,
  //   //   model_eval_t, scalar_t, mysizer, nonlin_solve_t>;
  //   // //      aux_res_pol_t, aux_jac_pol_t>;
  //   // aux_stepper_t auxStepperObj(appObj, nonls);//, auxResObj, auxJacObj);

  //   // now define the target stepper 
  //   // using res_pol_t = ode::policy::implicitBDF2IncrementResidual<
  //   //   state_t, residual_t, model_eval_t, scalar_t>;
  //   // res_pol_t resObj(y0);
  //   // using jac_pol_t = ode::policy::implicitBDF2IncrementJacobian<
  //   //   state_t, jac_t, model_eval_t, scalar_t>;
  //   // jac_pol_t jacObj(y0);

  //   // //stepper
  //   // using stepper_t = ode::implicitBDF2Stepper<
  //   //   state_t, residual_t, jac_t, scalar_t,
  //   //   model_eval_t, scalar_t, mysizer, nonlin_solve_t,
  //   //   aux_stepper_t>;//, res_pol_t, jac_pol_t>;
  //   // stepper_t stepperObj(appObj, nonls, auxStepperObj);//, resObj, jacObj);

  //   //stepper
  //   //using stepper_t = ode::implicitEulerStepper<
  //   using stepper_t = ode::implicitAdamsMoulton1Stepper<
  //     state_t, residual_t, jac_t, scalar_t,
  //     model_eval_t, scalar_t, nonlin_solve_t>;//, res_pol_t, jac_pol_t>;
  //   stepper_t stepperObj(appObj, nonls);//, resObj, jacObj);
    
  //   ode::integrateNSteps( stepperObj, y0, 0.0, dt, final_t/dt, collectorObj);
  //   // // // // // std::cout << collectorObj.getCount() << std::endl;
  //   std::cout << "Final solution " << std::endl;
  //   std::cout << "\n-----------" << std::endl;
  //   std::cout << *y0.data();
  //   // // // // for (int i=0; i<y0.size(); ++i)
  //   // // // //   std::cout << std::setprecision(14) << y0[i]  << " ";
  //   // // // // std::cout << std::endl;
  // }

  
    
  return 0;
}
