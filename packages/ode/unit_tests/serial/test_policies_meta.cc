
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

template<class T1, class T2, class T3, class T4, class T5>
struct testResRight :
  public ode::policy::explicitResidualPolicyBase<
  testResRight,T1,T2,T3,T4,T5>
{};


TEST(policies_meta, checkInheritance)
{
  using vecd = std::vector<double>;
  using veci = std::vector<int>;
  using model_t = model<int,double>;

  static_assert( !ode::meta::derivesFromExplicitResidualPolicyBase
  		 <testResWrong<vecd, vecd>>::value, "");

  static_assert( ode::meta::derivesFromExplicitResidualPolicyBase
   		 <testResRight<vecd,vecd,veci,double,int>>::value, "");
  
  static_assert( ode::meta::derivesFromExplicitResidualPolicyBase<
  		 ode::policy::explicitEulerStandardResidual<vecd, vecd,
		 model_t, double>
  		 >::value, "");

  static_assert( !ode::meta::isExplicitEulerResidualStandardPolicy
  		 <testResWrong<vecd, vecd>>::value, "");

  static_assert( !ode::meta::isExplicitEulerResidualStandardPolicy
   		 <testResRight<vecd,vecd,veci,double,int>>::value, "");
  
  static_assert( ode::meta::isExplicitEulerResidualStandardPolicy<
  		 ode::policy::explicitEulerStandardResidual<vecd, vecd,
		 model_t, double>
  		 >::value, "");
}
//--------------------------------------------------------------------
