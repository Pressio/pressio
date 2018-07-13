
#include "CORE_VECTOR"
#include "CORE_MATRIX"
#include "ODE_ALL"
#include "SOLVERS_EXP"
// app class
#include "apps_burgers1d_eigen.hpp"
#include "ode_observer.hpp"
//#include "../apps_helper_ode.hpp"

struct mysizer{
 using state_t = core::vector<apps::burgers1dEigen::state_type>;
 static size_t getSize(state_t & obj){
   return obj.size();
 };
  static void matchSize(const state_t & src,
			state_t & obj){
    obj.resize(src.size());
 };
};

template<typename T>
void printSol(std::string mess,
	      const T & y){
  std::cout << mess << std::endl;
  for (int i=0; i<y.size(); ++i)
    std::cout << std::setprecision(14) << y[i]  << " ";
  std::cout << std::endl;
}


int main(int argc, char *argv[])
{
  using native_state_t = apps::burgers1dEigen::state_type;
  using native_jac_t = apps::burgers1dEigen::jacobian_type;
  using scalar_t = apps::burgers1dEigen::scalar_type;
  using model_eval_t = apps::burgers1dEigen;

  // create app object
  Eigen::Vector3d mu(5.0, 0.02, 0.02);
  model_eval_t appObj(mu, 10);
  appObj.setup();

  // integrate in time startxbi5ng from y0
  scalar_t dt = 0.01;
  scalar_t final_t = 35;

  // wrap with core structures
  using state_t = core::vector<native_state_t>;
  using jac_t = core::matrix<native_jac_t>;
  using residual_t = state_t;
  native_state_t y0n = appObj.getInitialState();
  snapshot_collector collObj;

  state_t y0(y0n);

  // linear solver
  using lin_solve_t =
    solvers::experimental::linearSolver<jac_t, state_t, state_t>;
  lin_solve_t ls;
  // nonlinear solver
  using nonlin_solve_t =
    solvers::experimental::newtonRaphson<state_t, state_t,
  					 jac_t, lin_solve_t>;
  nonlin_solve_t nonls(ls);

  // stepper
  // using stepper_t = ode::implicitEulerStepper<
  //   state_t, residual_t, jac_t, scalar_t,
  //   model_eval_t, scalar_t, mysizer, nonlin_solve_t>;
  // stepper_t stepperObj(appObj, nonls);//, resObj, jacObj);
  
  // using stepper_t = ode::explicitEulerStepper<
  //   state_t, residual_t, scalar_t,
  //   model_eval_t, scalar_t, mysizer>;
  // stepper_t stepperObj(appObj);//, resObj, jacObj);

  // ode::integrateNSteps(stepperObj, y0, 0.0, dt, final_t/dt, collObj);
  // printSol("final", y0);


  ///////////////////////////

  
  // exstd_t::explicitStandardEuler(y0, 0.0, final_t, dt, appObj, collObj);
  // printSol("Exp Euler", y0);

  // *y0.data() = y0n;  
  // exstd_t::explicitStandardRK4(y0, 0.0, final_t, dt, appObj, collObj);
  // printSol("Exp RK4", y0);



  




  // using res_pol_t = ode::policy::incrementBasedResidual<
  //   state_t, residual_t, model_eval_t, scalar_t, mysizer>;
  // res_pol_t resObj(y0);

  // using jac_pol_t = ode::policy::incrementBasedJacobian<
  //   state_t, jac_t, model_eval_t, scalar_t, mysizer>;
  // jac_pol_t jacObj(y0);    
  
  // {
  //   *y0.data() = y0n;
  // //   //********************************************
  // //   // EXPLICIT EULER
  // //   //********************************************
  //   using stepper_t = ode::explicitEulerStepper<
  //     state_t, residual_t, scalar_t,
  //     model_eval_t, scalar_t, mysizer>;
  //   stepper_t stepperObj(appObj);
    
  //   ode::integrateNSteps(stepperObj, y0, 0.0, dt, final_t/dt, collObj);
  //   //    std::cout << collectorObj.getCount() << std::endl;
  //   std::cout << "Final solution " << std::endl;
  //   for (int i=0; i<y0.size(); ++i)
  //     std::cout << std::setprecision(14) << y0[i]  << " ";
  //   std::cout << std::endl;
  // }

  // {
  //   *y0.data() = y0n;
  //   //********************************************
  //   // EXPLICIT RK4
  //   //********************************************
  //   using stepper_t = ode::explicitRungeKutta4Stepper<
  //     state_t, residual_t, scalar_t,
  //     model_eval_t, scalar_t, mysizer>;
  //   stepper_t stepperObj(appObj);
    
  //   ode::integrateNSteps(stepperObj, y0, 0.0, dt, final_t/dt, collObj);
  //   //    std::cout << collectorObj.getCount() << std::endl;
  //   std::cout << "Final solution " << std::endl;
  //   for (int i=0; i<y0.size(); ++i)
  //     std::cout << std::setprecision(14) << y0[i]  << " ";
  //   std::cout << std::endl;
  // }
  
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

  //   // linear solver
  //   using lin_solve_t = solvers::experimental::linearSolver<jac_t,state_t,state_t>;
  //   lin_solve_t ls;    
  //   // nonlinear solver
  //   using nonlin_solve_t =
  //     solvers::experimental::newtonRaphson<state_t,state_t,jac_t,lin_solve_t>;
  //   nonlin_solve_t nonls(ls);

  //   //stepper
  //   using stepper_t = ode::implicitAdamsMoulton1Stepper<
  //     state_t, residual_t, jac_t, scalar_t,
  //     model_eval_t, scalar_t, mysizer, nonlin_solve_t>;//, res_pol_t, jac_pol_t>;
  //   stepper_t stepperObj(appObj, nonls);//, resObj, jacObj);
    
  //   ode::integrateNSteps( stepperObj, y0, 0.0, dt, final_t/dt, collectorObj);
  //   // // // // std::cout << collectorObj.getCount() << std::endl;
  //   // std::cout << "Final solution " << std::endl;
  //   std::cout << "\n-----------" << std::endl;
  //   std::cout << std::setprecision(14) << *y0.data();
  //   // // // for (int i=0; i<y0.size(); ++i)
  //   // // //   std::cout << std::setprecision(14) << y0[i]  << " ";
  //   // // // std::cout << std::endl;
  // }

  // {
  //   *y0.data() = y0n;

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
  //   using aux_stepper_t = ode::implicitEulerStepper<
  //     state_t, residual_t, jac_t, scalar_t,
  //     model_eval_t, scalar_t, mysizer, nonlin_solve_t>;
  //   //      aux_res_pol_t, aux_jac_pol_t>;
  //   aux_stepper_t auxStepperObj(appObj, nonls);//, auxResObj, auxJacObj);

  //   // now define the target stepper 
  //   // using res_pol_t = ode::policy::implicitBDF2IncrementResidual<
  //   //   state_t, residual_t, model_eval_t, scalar_t>;
  //   // res_pol_t resObj(y0);
  //   // using jac_pol_t = ode::policy::implicitBDF2IncrementJacobian<
  //   //   state_t, jac_t, model_eval_t, scalar_t>;
  //   // jac_pol_t jacObj(y0);

  //   // //stepper
  //   using stepper_t = ode::implicitBDF3Stepper<
  //     state_t, residual_t, jac_t, scalar_t,
  //     model_eval_t, scalar_t, mysizer, nonlin_solve_t,
  //     aux_stepper_t>;//, res_pol_t, jac_pol_t>;
  //   stepper_t stepperObj(appObj, nonls, auxStepperObj);//, resObj, jacObj);
    
  //   ode::integrateNSteps( stepperObj, y0, 0.0, dt, final_t/dt, collectorObj);
  //   // // // // // std::cout << collectorObj.getCount() << std::endl;
  //   std::cout << "Final solution " << std::endl;
  //   std::cout << "\n-----------" << std::endl;
  //   std::cout << std::setprecision(14) << *y0.data();
  //   // // // // for (int i=0; i<y0.size(); ++i)
  //   // // // //   std::cout << std::setprecision(14) << y0[i]  << " ";
  //   // // // std::cout << std::endl;
  // }

     
  return 0;
}
