
#include <gtest/gtest.h>
#include "pressio_ode.hpp"
#include "../reference_apps_for_testing.hpp"

struct MyApp{
  using scalar_type   = double;
  using state_type    = std::vector<scalar_type>;
  using velocity_type = state_type;

public:
  void velocity(const state_type & y,
		scalar_type t,
    velocity_type & R) const{
    R[0] = -10. * y[0];
    R[1] = -10. * y[1];
    R[2] = -10. * y[2];
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

  static void do_update(v_t & v, const double c,
  			const v_t & v0, const double a,
  			const v_t & v1, const double b){
    for (size_t i=0; i<v.size(); ++i)
      v[i] = c*v[i] + a*v0[i] + b*v1[i];
  }

  static void do_update(v_t & v,
  			const v_t & v0, const double a,
  			const v_t & v1, const double b){
    for (size_t i=0; i<v.size(); ++i)
      v[i] = a*v0[i] + b*v1[i];
  }

  static void do_update(v_t & v,
			const v_t & v1, const double b,
			const v_t & v2, const double c,
			const v_t & v3, const double d,
			const v_t & v4, const double e)
  {
    for (size_t i=0; i<v.size(); ++i)
      v[i] = b*v1[i] + c*v2[i] + d*v3[i] + e*v4[i];
  }

  static void do_update(v_t & v, const double a,
			const v_t & v1, const double b,
			const v_t & v2, const double c,
			const v_t & v3, const double d,
			const v_t & v4, const double e)
  {
    for (size_t i=0; i<v.size(); ++i)
      v[i] = a*v[i] + b*v1[i] + c*v2[i] + d*v3[i] + e*v4[i];
  }
};

struct myops{
  using update_op = updateOps;
};


TEST(ode_explicit_rk4, userDefinedOps){
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
    ode::ExplicitEnum::RungeKutta4, state_t, app_t, res_t,
    double, myops>;
  stepper_t stepperObj(y, appObj);

  // integrate in time
  double dt = 0.1;
  ode::integrateNSteps(stepperObj, y, 0.0, dt, 1);
  {
    std::cout << std::setprecision(14) 
    << (*y.data())[0] << " " 
    << (*y.data())[1] << " " 
    << (*y.data())[2] << std::endl;

    auto yptr2 = y.data();
    EXPECT_DOUBLE_EQ( (*yptr2)[0], 0.375);
    EXPECT_DOUBLE_EQ( (*yptr2)[1], 0.75);
    EXPECT_DOUBLE_EQ( (*yptr2)[2], 1.125);
  }
}
