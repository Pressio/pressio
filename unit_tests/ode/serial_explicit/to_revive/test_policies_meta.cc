
#include <gtest/gtest.h>
#include "policies/base/ode_explicit_policy_base.hpp"
#include "policies/ode_explicit_euler_standard_policy.hpp"
#include "policies/ode_explicit_policies_meta.hpp"


template<typename T1,typename T2>
struct model{
  void residual(const T1 &y, T2 &R){
  };
};


template<typename T1, typename T2>
struct testResWrong{
};

TEST(policies_meta, checkInheritance)
{
  using vecd = std::vector<double>;
  using veci = std::vector<int>;
  using model_t = model<int,double>;

  static_assert( !ode::meta::is_explicit_euler_velocity_standard_policy
  		 <testResWrong<vecd, vecd>>::value, "");
  
  static_assert( ode::meta::is_explicit_euler_velocity_standard_policy<
  		 ode::policy::explicitEulerStandardResidual<vecd, vecd,
		 model_t, double>
  		 >::value, "");
}
//--------------------------------------------------------------------
