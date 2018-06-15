#include <iomanip>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <type_traits>
#include <cmath>
#include <fstream>
#include <cassert>

// integrator
#include "integrators/ode_integrate_n_steps.hpp"
// euler stepper
#include "step_methods/ode_explicit_euler_stepper.hpp"
#include "policies/ode_explicit_euler_standard_policy.hpp"
// rk4 
#include "step_methods/ode_explicit_runge_kutta4_stepper.hpp"
#include "policies/ode_explicit_runge_kutta4_standard_policy.hpp"
// vector wrappers
// #include "vector/core_vector_serial_eigen.hpp"
// #include "vector/core_vector_serial_stdlib.hpp"
//#include "policies/base/ode_explicit_policy_base.hpp"


constexpr double sigma = 10.0;
constexpr double R = 28.0;
constexpr double b = 8.0 / 3.0;
using time_type = double;

struct stdvectorResizer{
  using vecD = std::vector<double>;
  // this has to be default constructible  
  void operator()(const vecD & source, vecD & dest){
    if ( dest.size()!=source.size() )
      dest.resize(source.size());
  };
};
//---------------------------------------------------------------------
//---------------------------------------------------------------------

struct lorenzStd{
  using state_t = std::vector<double>;
  using residual_t = state_t;
  void residual( const state_t &x , residual_t &dxdt , double t ){
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = -b * x[2] + x[0] * x[1];
  }
};
//---------------------------------------------------------------------
//---------------------------------------------------------------------

template<typename state_t>
struct coll{
  void operator()(size_t step, time_type t, const state_t & x){
    // if ( step % 100 == 0){
      std::cout << "step: " << step << " "
		<< "t = " << t << " " 
		<< x[0] << " "
		<< x[1] << " "
		<< x[2] << std::endl;
      //    }
  }
};
//---------------------------------------------------------------------
//---------------------------------------------------------------------


template<class state_t, class res_t,
	 class mod_t, class time_t, class T4>
class expEulerMine
  : public ode::policy::explicitResidualPolicyBase<
           expEulerMine,state_t,res_t,mod_t,time_t,T4>
{
public:
  expEulerMine(T4 perturb) : bump_(perturb){}
  expEulerMine() = delete;
  ~expEulerMine() = default;
  T4 bump_;
private:
  void computeImpl(const state_t & y, res_t & res,
		   mod_t & model, time_t t){
    model.residual(y, res, t);
    std::cout << "THIS IS MY ONW bump = " << t << std::endl;
  }  
private:
  friend ode::policy::explicitResidualPolicyBase<
         expEulerMine,state_t,res_t,mod_t,time_t,T4>;
};



int main(int argc, char *argv[])
{
  using state_t = lorenzStd::state_t;
  using rhs_t = lorenzStd::state_t;
  using scalar_t = double;
  using model_t = lorenzStd;
  
  time_type final_time = 0.10;
  time_type dt = 0.01;
  
  // create app obj
  model_t appObj;

  // stepper and integrator
  // using residual_policy_t = myOwnExEu<state_t,rhs_t,model_t,double>;
  // residual_policy_t resPol;

  {
    // initial condition
    state_t y0 = { 10.0 , 1.0 , 1.0 };

    // using stepper_t =
    //   ode::explicitEulerStepper<state_t, rhs_t, scalar_t,
    // 				stdvectorResizer, model_t,
    // 				time_type /*, default: standard policy residual */>;
    using stepper_t =
      ode::explicitRungeKutta4Stepper<state_t, rhs_t, scalar_t,
  				stdvectorResizer, model_t,
  				time_type /*, default: standard policy residual */>;
    stepper_t stepperObj(appObj);
    coll<state_t> observer;
    ode::integrateNSteps( stepperObj, y0, observer, 0.0, dt, final_time/dt );
  }
  // {
  //   // initial condition
  //   state_t y0 = { 10.0 , 1.0 , 1.0 }; 

  //   using residual_policy_t = expEulerMine<state_t, rhs_t,
  // 					   model_t, time_type, int>;

  //   static_assert(!ode::meta::isExplicitEulerResidualStandardPolicy<
  // 		  residual_policy_t>::value, "");
  //   using stepper_t =
  //     ode::explicitEulerStepper<state_t, rhs_t, scalar_t,stdvectorResizer,
  //   				model_t,time_type, residual_policy_t>;
  //   residual_policy_t resObj(55);
  //   stepper_t stepperObj(appObj, resObj);
  //   coll<state_t> observer;
  //   ode::integrateNSteps( stepperObj, y0, observer, 0.0, dt, final_time/dt );
  // }


  
  // {
  //   // omitting the policy defaults to creating a corresponding standard residual policy
  //   using stepper_t =
  //     ode::explicitRungeKutta4Stepper<state_t, rhs_t, scalar_t,
  // 				      stdvectorResizer, model_t /*, standard policy residual */>;
  //   stepper_t stepperObj(appObj);
  //   coll<state_t> observer;
  //   ode::integrateNSteps( stepperObj, y0, observer, 0.0, dt, final_time/dt );
  // }
  
  
  // do time integration wrapping with my vector
  //  using myvec_t = core::vector<state_t>;
  // myvec_t y(y0);
  //static_assert( core::details::traits<myvec_t>::isVector == 1, "ok");
  //static_assert( core::details::traits<myvec_t>::isEigen == 0, "ok"); 
  return 0;
}







// template<typename state_type, typename residual_type,
// 	 typename model_type, typename time_type>
// class myOwnExEu
//   : public ode::policy::explicitEulerResidualPolicyBase< myOwnExEu<state_type,
// 						      residual_type,
// 						      model_type,
// 						      time_type
// 						      >,
// 					    state_type, residual_type,
// 					    model_type, time_type
// 					    >
// {
//  public:
//   using derived_t = myOwnExEu<state_type, residual_type,
// 			      model_type, time_type >;
//   using base_t = ode::policy::explicitEulerResidualPolicyBase<derived_t,
//  						 state_type, residual_type,
//  						 model_type, time_type >;

//   void computeImpl(const state_type & y, residual_type & R,
//  		   model_type & model, time_type t){
//     std::cout << " evaluating " << std::endl;
//     model.residual(y, R, t);
//   }
// };






// using vecD = std::vector<double>;
// using ui_t = unsigned int;


// class burg1d
// {
// public:
//   using state_type = std::vector<double>;

// public:  
//   burg1d(vecD params) : mu_(params){}

//   void setup(){
//     dx_ = (xR_ - xL_)/static_cast<double>(Ncell_);
//     xGrid_.resize(Ncell_,0.0);
//     for (ui_t i=0; i<Ncell_; ++i){
//       xGrid_[i] = dx_*i + dx_*0.5;
//     };
//     U_.resize(Ncell_, 1.0);
//   };

//   state_type viewInitCond(){
//     return U_;
//   };
  
//   void operator() ( const state_type & u,
// 		    state_type & R,
// 		    const double /* t */ )
//   {
//     R[0] = (0.5 * ( mu_[0]*mu_[0] - u[0]*u[0] ) )/dx_;
//     for (ui_t i=1; i<Ncell_; ++i){
//       R[i] = ( 0.5*(u[i-1]*u[i-1] - u[i]*u[i]) )/dx_;
//     }
//     for (ui_t i=0; i<Ncell_; ++i){
//       R[i] += mu_[1]*exp(mu_[2]*xGrid_[i]);
//     }    
//   }

//   void integrate()
//   {
//     // int nsteps = 2500;
//     // vecD rhs(Ncell_,0.0);
//     // double dt = tfinal_/static_cast<double>(nsteps);

//     // std::ofstream file;
//     // file.open( "out.txt" );
//     // for (ui_t step=0; step<nsteps; ++step)
//     // {
//     //   (*this)(U_, rhs, step*dt);
//     //   for (ui_t i=0; i<Ncell_; ++i){
//     //   	U_[i] += dt*(rhs[i]);
//     // 	if (step % 50 == 0 || step==0)
//     // 	 file << std::fixed << std::setprecision(10) << U_[i] << " ";
//     //   }      
//     //   if (step % 50 == 0)
//     // 	file << std::endl;
//     // }
//     // file.close();
//   }

// private:
  
//   const double xL_ = 0.0;
//   const double xR_ = 100.0;
//   const ui_t Ncell_ = 10;
//   vecD mu_;
//   const double t0 = 0.0;
//   const double tfinal_ = 35.0;
//   double dx_;
//   vecD xGrid_;
//   state_type U_;
// };
// //==================================================

// struct stateResizer{
//   // has to be default constructible
  
//   void operator()(const vecD  & source,
// 		  vecD & dest)
//   {
//     if ( dest.size()!=source.size() )
//       dest.resize(source.size());
//   };
// };


// struct snapshot_collector{
//   using matrix = std::vector<vecD>;
//   using state_t = vecD;
  
//   matrix snaps_;
//   size_t count_;
//   void operator()(size_t step, double t, const state_t & x){
//     if (step % 50 ==0 ){
//       snaps_.emplace_back(x);
//       count_++;
//     }
//     std::cout << step << " " << t << std::endl;
//   }
//   size_t getCount(){
//     return count_;
//   };

//   void print()
//   {
//     std::ofstream file;
//     file.open( "out.txt" );
//     for (ui_t step=0; step<count_; ++step)
//     {
//       for (auto & it : snaps_[step])
//     	 file << std::fixed << std::setprecision(10) << it << " "; 
//       file << std::endl;
//     }
//     file.close();    
//   };
// };


// // template <typename state_type, typename oapp>
// // class GP()
// // {
// // private:
// //   state_type myY_;
// //   oapp * appPtr_;
    
// // public:
// //   GP(state_type & src, ...)
// //     : myY(src), ... {}
// //   ~GP(){}
  
// //   void operator()()
// //   {
// //     (*oappPtr)()( V' * y )
// //   }

// //   void run()
// //   {
// //     ode::eulerStepper<state_t,state_t,double,stateResizer> myStepper;
// //     ode::integrateNSteps(myStepper, *this, myY, snColl, 0.0, 0.0035, 10000);    
// //   }
// // };


// int main(int argc, char *argv[])
// {
//   const vecD mu{5.0, 0.02, 0.02};
//   burg1d app(mu);
//   app.setup();

//   using state_t = burg1d::state_type;
//   state_t U = app.viewInitCond();
//   snapshot_collector snColl;  
//   {
//     ode::eulerStepper<state_t,state_t,double,stateResizer> myStepper;
//     ode::integrateNSteps(myStepper, app, U, snColl, 0.0, 0.0035, 10000);
//     //ode::integrateNSteps(myStepper, app, U, 0.0, 0.0035, 10000);
//     // for (auto & it : U)
//     //   std::cout << it  << std::endl;
//     std::cout << snColl.getCount() << std::endl;
//     snColl.print();
//   }

//   // //----------------------------------------  
//   // // lets assume we have snapshot matrix V
  
//   // using myvec_t = core::vector<state_t>;
//   // myvec_t y(U); // y contains the initial condition
//   // y.matmul(V); // something like this

//   // GP gpSolver(y, app, ...);
//   // gpSolver.run();
  
//   // // auto const * data = gigi.view();
//   // // for (auto & it : *data)
//   // //   std::cout << it  << std::endl;
  
//   return 0;
// }
