
#include <gtest/gtest.h>
#include "pressio_ode.hpp"
#include "../reference_apps_for_testing.hpp"

namespace 
{

template<typename state_type>
class MyPolicyOk
{
public:
  template <typename system_type>
  state_type create(const system_type & system);

  template <typename system_type, typename scalar_type>
  void compute(const state_type & state,
	       state_type & f,
	       const system_type & system,
	       const scalar_type & time);
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
}

TEST(user_defined_model, admissibleVelocityPolicy)
{
  using namespace pressio;
  using app_t = ::pressio::ode::testing::refAppEigen;

  using state_w_type = pressio::containers::Vector<Eigen::VectorXd>;
  using policy_t = MyPolicyOk<state_w_type>;
  static_assert
  (ode::constraints::explicit_velocity_policy<policy_t, double, 
  	state_w_type, state_w_type, app_t>::value, "");
}

TEST(user_defined_model, nonAdmissibleVelocityPolicy)
{
  using namespace pressio;
  using app_t = ::pressio::ode::testing::refAppEigen;
  using state_w_type = pressio::containers::Vector<Eigen::VectorXd>;
  using policy_t = MyPolicyNotOk<state_w_type>;
  static_assert
  (!ode::constraints::explicit_velocity_policy<policy_t, double, 
  	state_w_type, state_w_type, app_t>::value, "");
}


TEST(user_defined_model, admissibleExplicitOde)
{
  using namespace pressio;
  using app_t = ::pressio::ode::testing::fakeAppForTraitsForExp;
  static_assert(
    ode::constraints::continuous_time_system_with_at_least_velocity<app_t>::value, "");
}
