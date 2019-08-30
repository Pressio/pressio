
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


struct updateOps{
  using v_t = std::vector<double>;

    static void do_update(v_t & v,
			  const v_t & v1, const double b){
    for (size_t i=0; i<v.size(); ++i)
      v[i] = b*v1[i];
  }

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
  auto yptr = y.data();
  EXPECT_DOUBLE_EQ( (*yptr)[0], 2.0);
  EXPECT_DOUBLE_EQ( (*yptr)[1], 4.);
  EXPECT_DOUBLE_EQ( (*yptr)[2], 6.);
  }
}
