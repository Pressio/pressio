
#include <gtest/gtest.h>
#include "pressio/ode_steppers_implicit.hpp"
#include "testing_apps.hpp"

namespace
{

struct ValidDiscreteTimeSystem{
  using scalar_type = float;
  using state_type = std::vector<float>;
  using discrete_time_residual_type = state_type;
  using discrete_time_jacobian_type = std::vector<std::vector<float>>;

  discrete_time_residual_type createDiscreteTimeResidual() const{ 
    return discrete_time_residual_type(); }
  discrete_time_jacobian_type createDiscreteTimeJacobian() const{ 
    return discrete_time_jacobian_type(); }

  template<class StepCountType>
  void discreteTimeResidual(StepCountType, 
                              double time, 
                              double dt, 
                              discrete_time_residual_type &, 
                              const state_type &) const{}

  template<class StepCountType>
  void discreteTimeResidual(StepCountType, 
                              double time, 
                              double dt, 
                              discrete_time_residual_type &, 
                              const state_type &,
                              const state_type &) const{}

  template<class StepCountType>
  void discreteTimeJacobian(StepCountType, 
                              double time, 
                              double dt, 
                              discrete_time_jacobian_type &, 
                              const state_type &) const{}

  template<class StepCountType>
  void discreteTimeJacobian(StepCountType, 
                              double time, 
                              double dt, 
                              discrete_time_jacobian_type &, 
                              const state_type &,
                              const state_type &) const{}
};

template<typename StateType, typename ResidualType>
class ResidualPolicy
{
public:
  using residual_type = ResidualType;

  residual_type create() const{ return residual_type(); }

  template <typename odetag, typename prev_states_type, class rhs_container>
  void compute(const StateType & y,
      const prev_states_type & oldYs,
      rhs_container & rhs,
      const double & t,
      const double & dt,
      int32_t step,
      residual_type & R) const
  {}
};//end class


template<typename StateType, typename JacobianType>
class JacobianPolicy
{
public:
  using jacobian_type = JacobianType;

  template <typename odetag, typename prev_states_type>
  void compute(const StateType & y,
      const prev_states_type & oldYs,
      const double &  t,
      const double &  dt,
      int32_t step,
      jacobian_type & J) const
  {}

  jacobian_type create() const{ return jacobian_type(); }
};
}

TEST(ode, concepts_discrete_time_system)
{
  using namespace pressio::ode;
  static_assert(discrete_time_system_with_user_provided_jacobian<ValidDiscreteTimeSystem>::value, "");
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
     residual_policy_t, state_t, double>::value, "");

  using jacobian_policy_t = JacobianPolicy<state_t, jac_t>;
  static_assert(ode::implicit_euler_jacobian_policy<
     jacobian_policy_t, state_t, double>::value, "");
}

TEST(ode, implicit_stencil_size)
{
  namespace po = pressio::ode;
  static_assert(po::implicit_stencil_size(po::implicitmethods::BDF1())==2,"" );
  static_assert(po::implicit_stencil_size(po::implicitmethods::BDF2())==3, "");
  static_assert(po::implicit_stencil_size(po::implicitmethods::CrankNicolson())==2,"");
}

