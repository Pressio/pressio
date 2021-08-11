
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

template<typename state_type>
class MyPolicyOk
{
public:
  template <typename system_type>
  state_type create(const system_type & system) const;

  template <typename system_type, typename scalar_type>
  void compute(const state_type & state,
	       state_type & f,
	       const system_type & system,
	       const scalar_type & time) const;
};

template<typename state_type>
class MyPolicyNotOk
{
public:
  // remove the create on purpose for test below
  // template <typename system_type>
  // state_type create(const system_type & system);

  template <typename system_type, typename scalar_type>
  void compute(const state_type & state,
	       state_type & f,
	       const system_type & system,
	       const scalar_type & time);
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

  using state_type = typename app_t::state_type;
  using time_type = double;
  static_assert(ode::explicit_velocity_policy<
    MyPolicyOk<state_type>, time_type, state_type, state_type, app_t>::value, "");

  static_assert(!ode::explicit_velocity_policy<
    MyPolicyNotOk<state_type>, time_type, state_type, state_type, app_t>::value, "");
}

TEST(ode, concepts_explicitly_steppable)
{
  using namespace pressio::ode;
  using state_type = std::vector<float>;
  static_assert(explicitly_steppable<SteppableClass, state_type, double>::value, "");
}
