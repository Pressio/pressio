
#include <gtest/gtest.h>

using ScalarType = double;
using VectorType = std::vector<ScalarType>;

namespace
{
struct MyApp{
  using independent_variable_type   = ScalarType;
  using state_type    = VectorType;
  using right_hand_side_type = VectorType;

public:
  state_type createState() const{ return state_type(3); }

  void rightHandSide(const state_type & y,
		independent_variable_type evalt,
                right_hand_side_type & R) const
  {
    R[0] = -10. * y[0];
    R[1] = -10. * y[1];
    R[2] = -10. * y[2];
  };

  right_hand_side_type createRightHandSide() const{
    right_hand_side_type R(3);
    return R;
  };
};
} // end namespace

#include "pressio/type_traits.hpp"

namespace pressio{

template<>
struct Traits<VectorType>
{
  using scalar_type = ScalarType;
};

namespace ops{
void deep_copy(VectorType & dest, const VectorType & from){
  dest = from;
}

VectorType clone(const VectorType & src)
{
  return VectorType(src.size());
}

void update(VectorType & v,        const ScalarType a,
		        const VectorType & v1, const ScalarType b)
{
  for (size_t i=0; i<v.size(); ++i){
    v[i] = a*v[i] + b*v1[i];
  }
}

void update(VectorType & v,        const ScalarType a,
            const VectorType & v1, const ScalarType b,
            const VectorType & v2, const ScalarType c)
{
  for (size_t i=0; i<v.size(); ++i){
    v[i] = a*v[i] + b*v1[i] + c*v2[i];
  }
}

void update(VectorType & v, const ScalarType a,
    const VectorType & v1, const ScalarType b,
    const VectorType & v2, const ScalarType c,
    const VectorType & v3, const ScalarType d,
    const VectorType & v4, const ScalarType e)
{
  for (size_t i=0; i<v.size(); ++i){
    v[i] = a*v[i] + b*v1[i] + c*v2[i] + d*v3[i] + e*v4[i];
  }
}
}} //end namespace pressio::ops

#include "pressio/ode_advancers.hpp"
#include "pressio/ode_steppers_explicit.hpp"

TEST(ode, explicit_euler_custom_types)
{
  using namespace pressio;
  using app_t	   = MyApp;
  using state_t = typename app_t::state_type;
  app_t appObj;

  state_t y(3);
  y[0] = 1.; y[1] = 2.; y[2] = 3.;
  auto stepperObj = ode::create_forward_euler_stepper(appObj);

  ScalarType dt = 0.1;
  ode::advance_n_steps(stepperObj, y, 0.0, dt, ::pressio::ode::StepCount(1));
  EXPECT_DOUBLE_EQ( y[0], 0.);
  EXPECT_DOUBLE_EQ( y[1], 0.);
  EXPECT_DOUBLE_EQ( y[2], 0.);
}

TEST(ode, explicit_rk4_custom_types)
{
  using namespace pressio;
  using app_t     = MyApp;
  using state_t    = typename app_t::state_type;
  app_t appObj;
  state_t y(3);
  y[0] = 1.; y[1] = 2.; y[2] = 3.;

  auto stepperObj = ode::create_runge_kutta4_stepper(appObj);

  ScalarType dt = 0.1;
  ode::advance_n_steps(stepperObj, y, 0.0, dt, ::pressio::ode::StepCount(1));
  EXPECT_DOUBLE_EQ( y[0], 0.375);
  EXPECT_DOUBLE_EQ( y[1], 0.75);
  EXPECT_DOUBLE_EQ( y[2], 1.125);
}


struct MyAppForAb2
{
  using independent_variable_type   = ScalarType;
  using state_type    = VectorType;
  using right_hand_side_type = VectorType;

public:
  state_type createState() const{
    return state_type(3);
  }

  void rightHandSide(const state_type & y,
		independent_variable_type evaltime,
		right_hand_side_type & R) const
  {
    R[0] = evaltime;
    R[1] = evaltime;
    R[2] = evaltime;
  };

  right_hand_side_type createRightHandSide() const
  {
    right_hand_side_type R(3);
    return R;
  };
};

struct CollectorTestAb2
{
  void operator()(const ::pressio::ode::StepCount & stepIn,
		  double time,
		  const VectorType & y)
  {
    const auto step = stepIn.get();
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

TEST(ode, explicit_ab2_custom_types)
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
  using app_t    = MyAppForAb2;
  using state_t = typename app_t::state_type;
  app_t appObj;
  state_t y(3);
  y[0] = 1.; y[1] = 2.; y[2] = 3.;

  auto stepperObj = ode::create_adams_bashforth2_stepper(appObj);
  ScalarType dt = 2.;
  CollectorTestAb2 C;
  ode::advance_n_steps(stepperObj, y, 0.0, dt, ::pressio::ode::StepCount(3), C);
}
