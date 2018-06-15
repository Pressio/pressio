
#include <gtest/gtest.h>
#include "./step_methods/ode_explicit_euler_stepper.hpp"


template<class T1, class T2, class T3, class T4>
struct fakePol :
    public ode::policy::explicitResidualPolicyBase<
      fakePol,T1,T2,T3,T4>
{
private:
  friend ode::policy::explicitResidualPolicyBase<
  fakePol,T1,T2,T3,T4>;
};


TEST(ode_explicit_euler_stepper, traits)
{
  using vecd = std::vector<double>;
  using state_t = vecd;
  using residual_t = vecd;
  using scalar_t = double;

  struct resizer{};
  struct app{};
  
  using resizer_t = resizer;
  using model_t = app;
  using time_type = double;
  using res_policy_t = fakePol<vecd,vecd,
			       model_t,time_type>;

  using stepper_t =
    ode::explicitEulerStepper<state_t, residual_t, scalar_t,
			      resizer_t, model_t, time_type,
			      res_policy_t>;

  static_assert(
  		!ode::meta::isExplicitEulerResidualStandardPolicy<
  		res_policy_t>::value,
  		"");
    
  app appObj;
  res_policy_t  polObj;
  stepper_t obj(appObj, polObj);

   using traits = ode::details::traits<stepper_t>;
  ::testing::StaticAssertTypeEq<typename
  				traits::state_t,vecd>();
  ::testing::StaticAssertTypeEq<typename
  				traits::residual_t,vecd>();
  ::testing::StaticAssertTypeEq<typename
  				traits::scalar_t,double>();
  ::testing::StaticAssertTypeEq<typename
  				traits::resizer_t,resizer>();
  ::testing::StaticAssertTypeEq<typename
  				traits::model_t,app>();
  ::testing::StaticAssertTypeEq<typename
  				traits::time_t,double>();
  ::testing::StaticAssertTypeEq<typename
  				traits::residual_policy_t,
  				res_policy_t>();
  static_assert( traits::order_value == 1, "");
}




TEST(ode_explicit_euler_stepper, traits2)
{
  using vecd = std::vector<double>;
  using state_t = vecd;
  using residual_t = vecd;
  using scalar_t = double;

  struct resizer{};
  struct app{};
 
  using resizer_t = resizer;
  using model_t = app;
  using time_type = double;

  using stepper_t =
    ode::explicitEulerStepper<state_t, residual_t, scalar_t,
			      resizer_t, model_t, time_type
			      /*res_policy defaulted*/>;
  app appObj;
  stepper_t obj(appObj);

  using traits = ode::details::traits<stepper_t>;
  ::testing::StaticAssertTypeEq<typename
  				traits::state_t,vecd>();
  ::testing::StaticAssertTypeEq<typename
  				traits::residual_t,vecd>();
  ::testing::StaticAssertTypeEq<typename
  				traits::scalar_t,double>();
  ::testing::StaticAssertTypeEq<typename
  				traits::resizer_t,resizer>();
  ::testing::StaticAssertTypeEq<typename
  				traits::model_t,app>();
  ::testing::StaticAssertTypeEq<typename
  				traits::time_t,double>();
  ::testing::StaticAssertTypeEq<typename
  				traits::residual_policy_t,
  				ode::policy::explicitEulerStandardResidual<
  				  state_t,residual_t,model_t,time_type>
  				>();
  static_assert( traits::order_value == 1, "");
}
