
#include <gtest/gtest.h>
#include "pressio/ode_implicit.hpp"
#include "testing_apps.hpp"

namespace
{

template<typename state_type, typename residual_type>
class ResidualPolicy
{
public:
  template<typename system_type>
  residual_type create(const system_type & model) const{ return residual_type(); }

  template <typename odetag, typename system_type, typename prev_states_type>
  void compute(const state_type & y,
      const prev_states_type & oldYs,
      const system_type & model,
      const double & t,
      const double & dt,
      ::pressio::ode::step_count_type step,
      residual_type & R) const
  {}
};//end class


template<typename state_type, typename jacobian_type>
class JacobianPolicy
{
public:
  template <typename odetag, typename system_type, typename prev_states_type>
  void compute(const state_type & y,
      const prev_states_type & oldYs,
      const system_type & model,
      const double &  t,
      const double &  dt,
      ::pressio::ode::step_count_type step,
      jacobian_type & J) const
  {}

  template<typename system_type>
  jacobian_type create(const system_type & model) const{
    return jacobian_type();
  }
};
}

TEST(ode, concepts_policies_arbitrary_stepper)
{
  using namespace pressio;

  using app_t = ode::testing::refAppForImpEigen;
  using state_t = typename app_t::state_type;
  using res_t = typename app_t::velocity_type;
  using jac_t = typename app_t::jacobian_type;

  using residual_policy_t = ResidualPolicy<state_t, res_t>;
  static_assert(ode::implicit_euler_residual_policy<
     residual_policy_t, state_t, res_t, app_t, double>::value, "");

  using jacobian_policy_t = JacobianPolicy<state_t, jac_t>;
  static_assert(ode::implicit_euler_jacobian_policy<
     jacobian_policy_t, state_t, jac_t, app_t, double>::value, "");
}

TEST(ode, implicit_stencil_size)
{
  namespace po = pressio::ode;
  static_assert(po::implicit_stencil_size(po::implicitmethods::BDF1())==2,"" );
  static_assert(po::implicit_stencil_size(po::implicitmethods::BDF2())==3, "");
  static_assert(po::implicit_stencil_size(po::implicitmethods::CrankNicolson())==2,"");
}

