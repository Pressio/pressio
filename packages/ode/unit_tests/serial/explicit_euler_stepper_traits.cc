
#include <gtest/gtest.h>
#include "CORE_ALL"
#include "ODE_ALL"

// template<class T1, class T2, class T3, class T4>
// struct fakePol :
//     public ode::policy::ExplicitResidualPolicyBase<
//       fakePol,T1,T2,T3,T4>
// {
// private:
//   friend ode::policy::ExplicitResidualPolicyBase<
//   fakePol,T1,T2,T3,T4>;
// };


TEST(ode_explicit_euler_stepper, traits)
{
  using namespace rompp;
  
  struct fakeapp{
    using scalar_type = double;
    using state_type = std::vector<double>;
    using space_residual_type = std::vector<double>;

    void residual(const state_type & y,
		  space_residual_type & R,
		  double t){
    };
    double residual(){
      return 0.0;
    };

  };
  fakeapp APP;

  using nstate_t = Eigen::VectorXd;
  using nresidual_t = Eigen::VectorXd;
  //  using scalar_t =double;
  using app_t = fakeapp;
  
  using state_t = core::Vector<nstate_t>;
  using res_t = core::Vector<nresidual_t>;
  state_t y0;
  res_t r0;
  
  using stepper_t = ode::ExplicitStepper<ode::ExplicitSteppersEnum::Euler,
  					 state_t, app_t>;
  stepper_t stepper(APP, y0, r0);

  // using stepper_t2 = ode::ExplicitStepper<ode::ExplicitSteppersEnum::Euler,
  // 					  app_t, state_t>;
  // stepper_t2 stepper2(APP, y0, r0);

  // using stepper_t3 = ode::ExplicitStepper<ode::ExplicitSteppersEnum::Euler,
  // 					  app_t, state_t, double>;
  // stepper_t3 stepper3(APP, y0, r0);
  

//   template<typename ode_state_type,
// 	 typename model_type,
// 	 typename ode_residual_type = ode_state_type,
// 	 typename enable = void>
// lass ExplicitEulerStepperStandardPolicy;

// template<typename ode_state_type,
// 	 typename residual_policy_type,
// 	 typename ode_residual_type = ode_state_type,
// 	 typename enable = void>
// class ExplicitEulerStepperArbitraryPolicy;

 
  // using vecd = std::<double>;
  // using state_t = vecd;
  // using residual_t = vecd;
  // using scalar_t = double;

  // struct app{};  
  // using model_t = app;
  // using time_type = double;
  // using res_policy_t = fakePol<vecd,vecd,
  // 			       model_t,time_type>;

  // using stepper_t =
  //   ode::ExplicitEulerStepper<state_t, residual_t, scalar_t,
  // 			      model_t, time_type, res_policy_t>;

  // static_assert(
  // 		!ode::meta::is_explicit_euler_residual_standard_policy<
  // 		res_policy_t>::value,
  // 		"");
    
  // app appObj;
  // res_policy_t  polObj;
  // stepper_t obj(appObj, polObj);

  //  using traits = ode::details::traits<stepper_t>;
  // ::testing::StaticAssertTypeEq<typename
  // 				traits::state_t,vecd>();
  // ::testing::StaticAssertTypeEq<typename
  // 				traits::residual_t,vecd>();
  // ::testing::StaticAssertTypeEq<typename
  // 				traits::scalar_t,double>();
  // ::testing::StaticAssertTypeEq<typename
  // 				traits::model_t,app>();
  // ::testing::StaticAssertTypeEq<typename
  // 				traits::time_t,double>();
  // ::testing::StaticAssertTypeEq<typename
  // 				traits::residual_policy_t,
  // 				res_policy_t>();
  // static_assert( traits::order_value == 1, "");
}




// TEST(ode_explicit_euler_stepper, traits2)
// {
//   using vecd = std::vector<double>;
//   using state_t = vecd;
//   using residual_t = vecd;
//   using scalar_t = double;

//   struct app{}; 
//   using model_t = app;
//   using time_type = double;

//   using stepper_t =
//     ode::ExplicitEulerStepper<state_t, residual_t, scalar_t,
// 			      model_t, time_type
// 			      /*res_policy defaulted*/>;
//   app appObj;
//   stepper_t obj(appObj);

//   using traits = ode::details::traits<stepper_t>;
//   ::testing::StaticAssertTypeEq<typename
//   				traits::state_t,vecd>();
//   ::testing::StaticAssertTypeEq<typename
//   				traits::residual_t,vecd>();
//   ::testing::StaticAssertTypeEq<typename
//   				traits::scalar_t,double>();
//   ::testing::StaticAssertTypeEq<typename
//   				traits::model_t,app>();
//   ::testing::StaticAssertTypeEq<typename
//   				traits::time_t,double>();
//   ::testing::StaticAssertTypeEq<typename
//   				traits::residual_policy_t,
//   				ode::policy::explicitEulerStandardResidual<
//   				  state_t,residual_t,model_t,time_type>
//   				>();
//   static_assert( traits::order_value == 1, "");
// }
