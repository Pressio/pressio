
#include <gtest/gtest.h>
#include "pressio/ode_steppers_explicit.hpp"

namespace{

struct ValidSystemWithVelocity1{
  using scalar_type = float;
  using state_type = std::vector<float>;
  using velocity_type = state_type;
  velocity_type createVelocity() const{ return velocity_type(); }
  void velocity(const state_type &, double time, velocity_type &) const{}
};

struct InvalidSystemWithVelocity1{
  using scalar_type = float;
  using state_type = std::vector<float>;
  using velocity_type = state_type;
  // velocity_type createVelocity() const{ return velocity_type(); }
  void velocity(const state_type &, double time, velocity_type &) const{}
};

struct InvalidSystemWithVelocity2{
  // using scalar_type = float;
  using state_type = std::vector<float>;
  using velocity_type = state_type;
  velocity_type createVelocity() const{ return velocity_type(); }
  void velocity(const state_type &, double time, velocity_type &) const{}
};

struct ValidSystemWithVelocityAndJacobian1{
  using scalar_type = float;
  using state_type = std::vector<float>;
  using velocity_type = state_type;
  using jacobian_type = std::vector<std::vector<float>>;
  velocity_type createVelocity() const{ return velocity_type(); }
  jacobian_type createJacobian() const{ return jacobian_type(); }
  void velocity(const state_type &, double time, velocity_type &) const{}
  void jacobian(const state_type &, double time, jacobian_type &) const{}
};
}//end anonym namespace 

TEST(ode, concepts_continuous_time_system)
{
  using namespace pressio::ode;
  static_assert(continuous_time_system_with_at_least_velocity<ValidSystemWithVelocity1>::value, "");
  static_assert(!continuous_time_system_with_at_least_velocity<InvalidSystemWithVelocity1>::value, "");
  static_assert(!continuous_time_system_with_at_least_velocity<InvalidSystemWithVelocity2>::value, "");
  static_assert(continuous_time_system_with_at_least_velocity<ValidSystemWithVelocityAndJacobian1>::value, "");

  static_assert(continuous_time_system_with_user_provided_jacobian<ValidSystemWithVelocityAndJacobian1>::value, "");
  static_assert(!continuous_time_system_with_user_provided_jacobian<ValidSystemWithVelocity1>::value, "");
}

