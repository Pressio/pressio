
#include <gtest/gtest.h>
#include "pressio/ode_advancers.hpp"

using ScalarType = double;
using VectorType = std::vector<ScalarType>;

struct MyFakeStepperExplicit
{
  void operator()(VectorType & odeState,
	      const ScalarType & time,
	      const ScalarType & dt,
	      const pressio::ode::step_count_type & step)  
  {
    for (std::size_t i=0; i<odeState.size(); i++){
      odeState[i] += time;
    }
  }
};

TEST(ode, explicit_advance_n_steps_fix_dt)
{
	VectorType odeState(5);
	std::for_each(odeState.begin(), odeState.end(), 
		[](ScalarType & val){val = 0.0; });
	MyFakeStepperExplicit stepper;
	ScalarType dt = 1.0;
	pressio::ode::advance_n_steps(stepper, odeState, 0.0, dt, 2.);
	std::for_each(odeState.begin(), odeState.end(), 
		[](const ScalarType & val){ EXPECT_DOUBLE_EQ(val, 1.);});
}


// template<typename ode_state_type>
// struct MyFakeStepper
// {
//   template<typename solver_type>
//   void operator()(ode_state_type & odeState,
// 	      const double & t,
// 	      const double & dt,
// 	      const pressio::ode::step_type & step,
// 	      solver_type & solver)
//   {
//     for (int i=0; i<odeState.extent(0); i++){
//       odeState(i) += dt;
//     }
//   }
// };

// struct MyFakeSolver
// {
//   template<typename system_t, typename state_t>
//   void solve(const system_t & sys, state_t & state){}
// };

// TEST(ode, advance_n_steps_with_dt_setter_and_collector)
// {
//   /*
//     at step 0: [1,2,3]
//     at step 1: dt=2,   and leads to [3,4,5]
//     at step 2: dt=0.5, and leads to [3.5,4.5,5.5]

//    */

//   using scalar_t = double;
//   using vec_t = Eigen::VectorXd;
//   using ode_state_t = pressio::containers::Vector<vec_t>;

//   ode_state_t y(3);
//   y(0) = 1.0; y(1) = 2.0; y(2) = 3.0;

//   MyFakeStepper<ode_state_t> stepper;
//   MyFakeSolver solver;

//   std::string checkStr= "PASSED";

//   auto dtManager = [](const ::pressio::ode::step_type & step,
// 		      const double & time,
// 		      double & dt)
// 		{
// 		  if(step==1) dt = 2.;
// 		  if(step==2) dt = 0.5;
// 		  if(step==3) dt = 1.;
// 		};

//   auto collector = [&checkStr](const ::pressio::ode::step_type & step,
// 		      const double & time,
// 		      const ode_state_t & y)
// 		   {
// 		     using vec_t = std::vector<scalar_t>;
// 		     vec_t true0 = {1., 2., 3.};
// 		     vec_t true1 = {3., 4., 5.};
// 		     vec_t true2 = {3.5, 4.5, 5.5};
// 		     vec_t true3 = {4.5, 5.5, 6.5};

// 		     vec_t * trueY = nullptr;

// 		     if (step==0) trueY = &true0;
// 		     if (step==1) trueY = &true1;
// 		     if (step==2) trueY = &true2;
// 		     if (step==3) trueY = &true3;

// 		     for (auto i=0; i<3; ++i){
// 		       if( std::abs(y(i)-(*trueY)[i]) > 1e-11 )
// 			 checkStr = "FAILED";
// 		     }
// 		   };
//   pressio::ode::advance_n_steps(stepper, y, 0., 3, dtManager, collector, solver);
// }
