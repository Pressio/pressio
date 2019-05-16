
#include <gtest/gtest.h>
#include "CORE_ALL"
#include "ODE_ALL"
#include "reference_apps_for_testing.hpp"


struct MyApp{

  using scalar_type   = double;
  using state_type    = std::vector<scalar_type>;
  using residual_type = state_type;

public:
  void residual(const state_type & y,
		residual_type & R,
		scalar_type t) const{
    R[0] = 10. * y[0];
    R[1] = 10. * y[1];
    R[2] = 10. * y[2];
  };

  residual_type residual(const state_type & y,
			 scalar_type t) const{
    residual_type R(y);
    residual(y, R, t);
    return R;
  };
};


struct updateOps{
  using v_t = std::vector<double>;

  static void do_update(v_t & v, const double a,
			const v_t & v1, const double b){
    for (size_t i=0; i<v.size(); ++i)
      v[i] = a*v[i] + b*v1[i];
  }
};

struct myops{
  using update_op = updateOps;
};


TEST(ode_explicit_euler, userDefinedOps){
  using namespace rompp;
  using app_t	    = MyApp;
  using nstate_t    = typename app_t::state_type;
  using nresidual_t = typename app_t::residual_type;
  app_t appObj;

  using state_t = core::Vector<nstate_t>;
  using res_t = core::Vector<nresidual_t>;
  state_t y(3);

  auto yptr = y.data();
  (*yptr)[0] = 1.; (*yptr)[1] = 2.; (*yptr)[2] = 3.;
  res_t r(3);
  appObj.residual(*y.data(), *r.data(), 0.0);

  using stepper_t = ode::ExplicitStepper<
    ode::ExplicitEnum::Euler, state_t, app_t, res_t,
    double, myops>;
  stepper_t stepperObj(y, appObj, r);

  // integrate in time
  double dt = 0.1;
  ode::integrateNSteps(stepperObj, y, 0.0, dt, 1ul);
  {
  auto yptr = y.data();
  EXPECT_DOUBLE_EQ( (*yptr)[0], 2.0);
  EXPECT_DOUBLE_EQ( (*yptr)[1], 4.);
  EXPECT_DOUBLE_EQ( (*yptr)[2], 6.);
  }
}
