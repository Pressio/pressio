
#include <gtest/gtest.h>
#include "policies/base/ode_explicit_policy_base.hpp"
#include "policies/ode_explicit_euler_standard_policy.hpp"
#include "policies/ode_explicit_policies_meta.hpp"


template<typename stt, typename rest>
struct app{
  void residual(const stt &y, rest &R, ode::details::time_type t){
    if (R.size()!=y.size())
      R.resize(y.size());
      
    for (size_t i=0; i<y.size(); i++)
      R[i] = y[i]*2 + t;
  };
};

TEST(explicit_euler_standard_policy, simpleTest)
{
  using state_type = std::vector<double>;
  using residual_type = std::vector<double>;
  using model_type = app<state_type, residual_type>;
  
  using pol_t =
    ode::policy::explicitEulerStandardResidual<state_type,residual_type,
					       model_type,double>;

  state_type y {2,4.};
  residual_type R;
  model_type appO;
    
  pol_t policyObj;
  policyObj.compute(y,R,appO,1);
  EXPECT_DOUBLE_EQ( R[0], 5.);
  EXPECT_DOUBLE_EQ( R[1], 9.);  
}
//--------------------------------------------------------------------
//--------------------------------------------------------------------


// flag_t is used here just for testing that we can have a policy that
// takes any # of temp arguments, since we want this to be as flexible as possible
template<class state_t, class res_t, class mod_t, class time_t, class flag_t>
class expEulerFakeNewPol
  : public ode::policy::ExplicitVelocityPolicyBase<expEulerFakeNewPol,state_t,
   						   res_t,mod_t,time_t,flag_t>{
public:
  expEulerFakeNewPol(flag_t perturb) : perturbation_(perturb){}
  ~expEulerFakeNewPol() = default;
  flag_t perturbation_;
private:
  void computeImpl(const state_t & y, res_t & R, mod_t & model, double t){    
    model.residual(y, R, t);
    R[0] += perturbation_[0];
    R[1] += perturbation_[1];
  }  
private:
  friend ode::policy::ExplicitVelocityPolicyBase<expEulerFakeNewPol,state_t,
						 res_t,mod_t,time_t,flag_t>;
};


TEST(explicit_euler_standard_policy, modifyStandard)
{
  using state_type = std::vector<double>;
  using residual_type = std::vector<double>;
  using model_type = app<state_type, residual_type>;
  using flag_t = std::map<int, double>;

  state_type y {2.,4.,1.};
  residual_type R;
  model_type appO;
  flag_t pp;
  pp[0] = -100.;
  pp[1] = -200.;

  using pol_t = expEulerFakeNewPol<state_type, residual_type,
				   model_type, double, flag_t>;
  pol_t policyObj(pp);
  policyObj.compute(y,R,appO,1);
  EXPECT_DOUBLE_EQ( R[0], -95.);
  EXPECT_DOUBLE_EQ( R[1], -191.);
}

