
#include <gtest/gtest.h>
#include "CORE_ALL"
#include "ODE_ALL"
#include "reference_apps_for_testing.hpp"


TEST(ode_explicit_euler, traits){
  using namespace rompp;
  
  using app_t = ode::testing::fakeAppForTraits;
  using nstate_t = typename app_t::state_type;
  using nres_t = typename app_t::residual_type;
  using state_t = core::Vector<nstate_t>;
  using res_t = core::Vector<nres_t>;
  
  using stepper_t = ode::ExplicitStepper<
    ode::ExplicitSteppersEnum::Euler, state_t, res_t, app_t>;

  using impl_t = typename stepper_t::base_t;
  using traits = ode::details::traits<impl_t>;
  
  static_assert(ode::meta::is_explicit_euler_residual_standard_policy<
  		typename traits::residual_policy_t>::value,
  		"");
  
  ::testing::StaticAssertTypeEq<typename
  				traits::state_t, state_t>();
  ::testing::StaticAssertTypeEq<typename
  				traits::residual_t,res_t>();
  ::testing::StaticAssertTypeEq<typename
  				traits::scalar_t,double>();
  ::testing::StaticAssertTypeEq<typename
  				traits::model_t,app_t>();
  static_assert( traits::order_value == 1, "");
}


TEST(ode_explicit_euler, numerics){
  using namespace rompp;
  
  using app_t = ode::testing::refAppEigen;
  using nstate_t = typename app_t::state_type;
  using nresidual_t = typename app_t::residual_type;
  app_t appObj;
  
  using state_t = core::Vector<nstate_t>;
  using res_t = core::Vector<nresidual_t>;
  state_t y(3);
  y[0] = 1.; y[1] = 2.; y[2] = 3.;
  res_t r(3);
  appObj.residual(*y.data(), *r.data(), 0.0);
  
  using stepper_t = ode::ExplicitStepper<
    ode::ExplicitSteppersEnum::Euler, state_t, app_t>;
  stepper_t stepperObj(appObj, y, r);
    
  // integrate in time 
  double dt = 0.1;
  ode::integrateNSteps(stepperObj, y, 0.0, dt, 1ul);
  EXPECT_DOUBLE_EQ( y[0], 1.1);
  EXPECT_DOUBLE_EQ( y[1], 2.2);
  EXPECT_DOUBLE_EQ( y[2], 3.3);
  //  std::cout << std::setprecision(14) << *y.data();

  // integrate in time 
  ode::integrateNSteps(stepperObj, y, 0.0, dt, 1ul);
  EXPECT_DOUBLE_EQ( y[0], 1.21);
  EXPECT_DOUBLE_EQ( y[1], 2.42);
  EXPECT_DOUBLE_EQ( y[2], 3.63);
}



