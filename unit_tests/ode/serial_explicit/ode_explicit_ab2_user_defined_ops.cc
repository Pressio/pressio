
#include <gtest/gtest.h>
#include "pressio_ode.hpp"

struct MyApp{
  using scalar_type   = double;
  using state_type    = std::vector<scalar_type>;
  using velocity_type = state_type;

public:
  void velocity(const state_type & y,
		scalar_type t,
		velocity_type & R) const
  {
    R[0] = t;
    R[1] = t;
    R[2] = t;
  };

  velocity_type createVelocity() const
  {
    velocity_type R(3);
    return R;
  };
};


struct Collector
{
  void operator()(const ::pressio::ode::types::step_t & step,
		  const double & time,
		  const std::vector<double> & y)
  {
    if (step==0){
      EXPECT_DOUBLE_EQ( y[0], 1.);
      EXPECT_DOUBLE_EQ( y[1], 2.);
      EXPECT_DOUBLE_EQ( y[2], 3.);
    }

    if (step==1){
      EXPECT_DOUBLE_EQ( y[0], 1.);
      EXPECT_DOUBLE_EQ( y[1], 2.);
      EXPECT_DOUBLE_EQ( y[2], 3.);
    }

    if (step==2){
      EXPECT_DOUBLE_EQ( y[0], 7.);
      EXPECT_DOUBLE_EQ( y[1], 8.);
      EXPECT_DOUBLE_EQ( y[2], 9.);
    }

    if (step==3){
      EXPECT_DOUBLE_EQ( y[0], 17.);
      EXPECT_DOUBLE_EQ( y[1], 18.);
      EXPECT_DOUBLE_EQ( y[2], 19.);
    }
  }
};

template <typename scalar_t>
struct MyOps
{
  using v_t = std::vector<scalar_t>;

  void deep_copy(v_t & dest, const v_t & from) const
  {
    dest = from;
  }

  void update(v_t & v, const v_t & v1, const scalar_t b) const
  {
    for (size_t i=0; i<v.size(); ++i)
      v[i] = b*v1[i];
  }

  void update(v_t & v, const scalar_t a,
	      const v_t & v1, const scalar_t b) const
  {
    for (size_t i=0; i<v.size(); ++i)
      v[i] = a*v[i] + b*v1[i];
  }

  void update(v_t & v,
	      const v_t & v1, const scalar_t b,
	      const v_t & v2, const scalar_t c) const
  {
    for (size_t i=0; i<v.size(); ++i)
      v[i] = b*v1[i] + c*v2[i];
  }

  void update(v_t & v, const scalar_t a,
	      const v_t & v1, const scalar_t b,
	      const v_t & v2, const scalar_t c) const
  {
    for (size_t i=0; i<v.size(); ++i)
      v[i] = a*v[i] + b*v1[i] + c*v2[i];
  }
};

TEST(ode_explicit_ab2, userDefinedOps)
{
  /*
    dy/dt = f
    where f returns [t t t]

    for AB2, step1 from t_0 -> t_1 is euler:
    y_1 = [1 2 3] + 2.*f(t0)
    f0 = [0 0 0]

    step2: from t_1 -> t_2
    y_2 = y_1 + dt*[ (3/2)*f(t_1) - (1/2)*f(t_0)]
        = [1 2 3] + 2*(3/2)*[2 2 2]
	= [7 8 9]

    step3:
    y_3 = y_2 + dt*[ (3/2)*f(t_2) - (1/2)*f(t_1)]
        = [7 8 9] + 2*[ (3/2)*[4 4 4] - (1/2)[2 2 2] ]
	= [17 18 19]
   */

  using namespace pressio;
  using app_t	   = MyApp;
  using scalar_t = typename app_t::scalar_type;
  using nstate_t = typename app_t::state_type;
  app_t appObj;

  using state_t = containers::Vector<nstate_t>;
  state_t y(3);

  auto yptr = y.data();
  (*yptr)[0] = 1.; (*yptr)[1] = 2.; (*yptr)[2] = 3.;

  using ops_t = MyOps<scalar_t>;
  ops_t myOps;

  auto stepperObj = ode::createAdamsBashforth2Stepper(y, appObj, myOps);
  double dt = 2.;
  Collector C;
  ode::advanceNSteps(stepperObj, y, 0.0, dt, 3, C);
}
