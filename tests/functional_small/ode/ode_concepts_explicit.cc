
#include <gtest/gtest.h>
#include "pressio/ode_explicit.hpp"

namespace
{

struct ValidSystemWithVelocity{
  using scalar_type = float;
  using state_type = std::vector<float>;
  using velocity_type = state_type;
  velocity_type createVelocity() const{ return velocity_type(); }
  void velocity(const state_type &, double time, velocity_type &) const{}
};

struct SteppableClass
{
  void doStep(const std::vector<float> state,
              const double time, 
              const double dt, 
              const pressio::ode::step_count_type step)
  {}
};

}//end anom namespace 

TEST(ode, concepts_velocityPolicy)
{
  using namespace pressio;
  using app_t = ValidSystemWithVelocity;

  static_assert
  (::pressio::ode::continuous_time_system_with_at_least_velocity<ValidSystemWithVelocity>::value,
   "");
}

TEST(ode, concepts_explicitly_steppable)
{
  using namespace pressio::ode;
  using state_type = std::vector<float>;
  static_assert(explicitly_steppable<SteppableClass, state_type, double>::value, "");
}
