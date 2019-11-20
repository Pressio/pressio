// this is the same as in other file
#include "wls_policies.hpp"
namespace pressio{ namespace rom{ namespace wls{
template<typename fom_type, typename wls_state_type, typename jtj_jtr_pol_t>
class WlsSystemJTJApi{
private:
  const jtj_jtr_pol_t & jtj_jtr_polObj_;
//  const fom_type & fomObj_;
public:
  using jtj_type = typename jtj_jtr_pol_t::jtj_t;
  using jtr_type = typename jtj_jtr_pol_t::jtr_t;
  WlsSystemJTJApi(const jtj_jtr_pol_t & jtj_jtr_polObj) : jtj_jtr_polObj_(jtj_jtr_polObj){}
//
  void JTJ_JTR(wls_state_type & wls_state, jtj_type & jtj, jtr_type & jtr) const{
//    // call policy, e.g. they have an evaluate method
//    //polObj_.evaluate( /*whatever goes in here */ )
    jtj_jtr_polObj_(wls_state,jtj,jtr);
}
};


// this is the same as in other file
template<typename fom_type, typename wls_state_type, typename res_pol_t, typename jac_pol_t>
class WlsSystemDefaultApi{
private:
  const res_pol_t & resPolObj_;
  const jac_pol_t & jacPolObj_;
public:
  using residual_type = typename res_pol_t::residual_t;
  using jacobian_type = typename jac_pol_t::jacobian_t;

  using fom_native_state_t      = typename fom_type::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using fom_container_type = typename pressio::ode::StatesContainer<fom_state_t, 2>;

  WlsSystemDefaultApi(const res_pol_t & resPolObj, const jac_pol_t & jacPolObj, int numSteps)
    : resPolObj_(resPolObj), jacPolObj_(jacPolObj), {}

  void residual(wls_state_type & wls_state, residual_type & resid) const{
      resPolObj_(wls_state,resid);
  }

  void jacobian(const wls_state_type & yROM, jacobian_type & Jphi) const{
//    jacobian_policy()(wls_state,fomObj_);
  };
};
}}}





