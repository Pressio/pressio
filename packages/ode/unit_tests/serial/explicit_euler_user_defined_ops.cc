
#include <gtest/gtest.h>
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "reference_apps_for_testing.hpp"


struct MyApp{
  using scalar_type   = double;
  using state_type    = std::vector<scalar_type>;
  using velocity_type = state_type;

public:
  void velocity(const state_type & y,
		scalar_type t, velocity_type & R) const{
    R[0] = 10. * y[0];
    R[1] = 10. * y[1];
    R[2] = 10. * y[2];
  };

  velocity_type velocity(const state_type & y,
			 scalar_type t) const{
    velocity_type R(y);
    velocity(y, t, R);
    return R;
  };
};

<template scalar_t>
struct updateOps{
  using v_t = std::vector<scalar_t>;

  static void do_update(v_t & v, const v_t & v1, const scalar_t b){
    for (size_t i=0; i<v.size(); ++i)
      v[i] = b*v1[i];
  }

  static void do_update(v_t & v, const scalar_t a,
			const v_t & v1, const scalar_t b){
    for (size_t i=0; i<v.size(); ++i)
      v[i] = a*v[i] + b*v1[i];
  }
};

<template scalar_t>
struct myops{
  // update_op is all you need to provide for explicit
  // time integration. This type that pressio will detect for doing
  // operations like vector additions.
  using update_op = updateOps<scalar_r>;

  // ... this might contains other types defining how to do
  // other operations different in nature.
  // for example how to do mat-vec products.
  // this will be addressed in a later tutorial.
};


TEST(ode_explicit_euler, userDefinedOps){
  using namespace pressio;
  using app_t	    = MyApp;
  using nstate_t    = typename app_t::state_type;
  using nveloc_t = typename app_t::velocity_type;
  app_t appObj;

  using state_t = containers::Vector<nstate_t>;
  using res_t = containers::Vector<nveloc_t>;
  state_t y(3);

  auto yptr = y.data();
  (*yptr)[0] = 1.; (*yptr)[1] = 2.; (*yptr)[2] = 3.;

  using stepper_t = ode::ExplicitStepper<
    ode::ExplicitEnum::Euler, state_t, app_t, res_t,
    double, myops>;
  stepper_t stepperObj(y, appObj);

  // integrate in time
  double dt = 0.1;
  ode::integrateNSteps(stepperObj, y, 0.0, dt, 1ul);
  {
  auto yptr2 = y.data();
  EXPECT_DOUBLE_EQ( (*yptr2)[0], 2.0);
  EXPECT_DOUBLE_EQ( (*yptr2)[1], 4.);
  EXPECT_DOUBLE_EQ( (*yptr2)[2], 6.);
  }
}
