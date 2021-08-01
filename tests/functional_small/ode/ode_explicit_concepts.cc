
#include <gtest/gtest.h>
#include "pressio_ode_explicit.hpp"
#include "testing_apps.hpp"

namespace
{

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
}

TEST(ode, admissibleVelocityPolicy)
{
  using namespace pressio;
  using app_t = ::pressio::ode::testing::refAppEigen;

  using state_type = Eigen::VectorXd;
  using policy_t = MyPolicyOk<state_type>;
  static_assert
  (pressio::ode::explicit_velocity_policy<policy_t, double,	state_type, state_type, app_t>::value, "");
}

TEST(ode, nonAdmissibleVelocityPolicy)
{
  using namespace pressio;
  using app_t = ::pressio::ode::testing::refAppEigen;
  using state_type = Eigen::VectorXd;
  using policy_t = MyPolicyNotOk<state_type>;
  static_assert
  (!ode::explicit_velocity_policy<policy_t, double, state_type, state_type, app_t>::value, "");
}

TEST(ode, admissibleExplicitOde)
{
  using namespace pressio;
  using app_t = ::pressio::ode::testing::fakeAppForTraitsForExp;
  static_assert(
    ode::continuous_time_system_with_at_least_velocity<app_t>::value, "");
}
