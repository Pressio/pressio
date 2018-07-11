
// #include "ode_observer.hpp"
// app class
#include "apps_burgers1d_epetra.hpp"
// integrator
#include "integrators/ode_integrate_n_steps.hpp"
// steppers
#include "explicit_steppers/ode_explicit_euler_stepper.hpp"
//#include "step_methods/ode_explicit_runge_kutta4_stepper.hpp"
#include "implicit_steppers/ode_implicit_euler_stepper.hpp"


// template<class state_t, class res_t,
// 	 class mod_t, class time_t, class T4>
// class expEulerMine
//   : public ode::policy::explicitResidualPolicyBase<
//            expEulerMine,state_t,res_t,mod_t,time_t,T4>
// {
// public:
//   expEulerMine(T4 perturb) : bump_(perturb){}
//   expEulerMine() = delete;
//   ~expEulerMine() = default;
//   T4 bump_;
// private:
//   void computeImpl(const state_t & y, res_t & res,
// 		   mod_t & model, time_t t){
//     model.residual(y, res, t);
//     std::cout << "THIS IS MY ONW bump = " << t << std::endl;
//   }  
// private:
//   friend ode::policy::explicitResidualPolicyBase<
//   expEulerMine,state_t,res_t,mod_t,time_t,T4>;
// };


struct mysizer{
  using state_t = apps::burgers1dEpetra::state_type;
 static size_t getSize(const state_t & obj){
   return obj.MyLength();
 };
 static void resize(state_t & obj, size_t newSize){};
};

int main(int argc, char *argv[])
{
  using state_t = apps::burgers1dEpetra::state_type;
  using jac_t = apps::burgers1dEpetra::jacobian_type;
  using scalar_t = apps::burgers1dEpetra::scalar_type;
  using target_app_t = apps::burgers1dEpetra;

  MPI_Init(&argc,&argv);
  int rank; // My process ID
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

  std::vector<double> mu({5.0, 0.02, 0.02});
  target_app_t appObj(mu, &Comm, 20);
  appObj.setup();
  
  // integrate in time starting from y0
  scalar_t final_time = 35.;
  scalar_t dt = 0.001;
  
  {
    //********************************************
    // EXPLICIT EULER
    //********************************************
    state_t y0 = appObj.getInitialState();

    auto const & dmap = y0.Map();
    using stepper_t = ode::explicitEulerStepper<
      state_t, state_t, scalar_t, target_app_t,
      scalar_t, mysizer>;
    stepper_t stepperObj(appObj, dmap);

    // using newPol_t = expEulerMine<state_t, state_t,
    // 				  target_app_t,
    // 				  double, double>;
    // newPol_t polObj(44);
      
    // using stepper_t2 = ode::explicitEulerStepper<
    //   state_t, state_t, scalar_t, target_app_t, scalar_t, newPol_t>;
    // stepper_t2 stepperObj2(appObj, polObj, dmap);
    
    ode::integrateNSteps(stepperObj, y0, 0.0, dt, final_time/dt);

    y0.Print(std::cout);
    
    //final_time/dt );
    // std::cout << collectorObj.getCount() << std::endl;
    // //    collectorObj.print();

    // std::cout << "Final solution " << std::endl;
    // for (int i=0; i<y0.size(); ++i)
    //   std::cout << std::setprecision(14) << y0[i]  << " ";
    // std::cout << std::endl;
  }
  
  MPI_Finalize();
  return 0;
}
