// this is the same as in other file
#include "wls_policies.hpp"
namespace pressio{ namespace rom{ namespace wls{
template<typename fom_type, typename wls_state_type, typename decoder_t, typename fom_state_container_t, typename jtj_jtr_pol_t>
class WlsSystemJTJApi{
private:
  const jtj_jtr_pol_t & jtj_jtr_polObj_;
//  const fom_type & fomObj_;
public:
  using scalar_t	= typename fom_type::scalar_type;
  using fom_native_state_t      = typename fom_type::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;

  fom_type & appObj_;
  

  using jtj_type = typename jtj_jtr_pol_t::hess_t;
  using jtr_type = typename jtj_jtr_pol_t::jtr_t;
  using fom_state_reconstr_t    = pressio::rom::FomStateReconstructor<fom_state_t, decoder_t>;
  using decoder_jac_t = typename decoder_t::jacobian_t;
  

  fom_state_reconstr_t  fomStateReconstructorObj_;
  fom_state_reconstr_t & fomStateReconstructorObjPtr;
  fom_state_container_t  fomStateContainerObj_;
  fom_state_container_t & fomStateContainerObjPtr;
  scalar_t dt_;
  std::size_t numStepsInWindow_; 
  WlsSystemJTJApi(fom_type & appObj, const jtj_jtr_pol_t & jtj_jtr_polObj, const decoder_t & decoderObj, fom_state_t & yFOM, scalar_t dt, std::size_t numStepsInWindow) : appObj_(appObj), jtj_jtr_polObj_(jtj_jtr_polObj),fomStateReconstructorObj_(yFOM,decoderObj),fomStateReconstructorObjPtr(fomStateReconstructorObj_),fomStateContainerObj_(yFOM),fomStateContainerObjPtr(fomStateContainerObj_)
{
  dt_ = dt;
  numStepsInWindow_ = numStepsInWindow;
}

  void JTJ_JTR(wls_state_type & wls_state, wls_state_type & wls_state_ic, jtj_type & jtj, jtr_type & jtr) const{
    jtj_jtr_polObj_(appObj_,wls_state,wls_state_ic,jtj,jtr,fomStateContainerObjPtr,fomStateReconstructorObjPtr,dt_,numStepsInWindow_);
}
};







// this is the same as in other file
template<typename fom_type, typename wls_state_type, typename res_pol_t, typename jac_pol_t, typename decoder_t>
class WlsSystemDefaultApi{
private:
  const res_pol_t & resPolObj_;
  const jac_pol_t & jacPolObj_;

public:
  using residual_type = typename res_pol_t::residual_t;
  using jacobian_type = typename jac_pol_t::jacobian_t;
  using scalar_t	= typename fom_type::scalar_type;

  using fom_native_state_t      = typename fom_type::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using fom_state_container_t = typename pressio::ode::StatesContainer<fom_state_t, 2>;
  using fom_state_reconstr_t    = pressio::rom::FomStateReconstructor<fom_state_t, decoder_t>;
  fom_state_reconstr_t  fomStateReconstructorObj_;
  fom_state_reconstr_t & fomStateReconstructorObjPtr;

  fom_state_container_t  fomStateContainerObj_;
  fom_state_container_t & fomStateContainerObjPtr;
  fom_type & appObj_;
  scalar_t dt_;

  WlsSystemDefaultApi(fom_type & appObj, const res_pol_t & resPolObj, const jac_pol_t & jacPolObj, const decoder_t & decoderObj, fom_state_t & yFOM, scalar_t dt)
    : appObj_(appObj), resPolObj_(resPolObj), jacPolObj_(jacPolObj), fomStateReconstructorObj_(yFOM,decoderObj),fomStateContainerObj_(yFOM), fomStateReconstructorObjPtr(fomStateReconstructorObj_),fomStateContainerObjPtr(fomStateContainerObj_)
  {
    dt_ = dt;
  }

  void residual(wls_state_type & wls_state, residual_type & resid) const{
      resPolObj_(appObj_,wls_state,resid,fomStateContainerObjPtr,fomStateReconstructorObjPtr,dt_);
  }

  void jacobian(const wls_state_type & yROM, jacobian_type & Jphi) const{
//    jacobian_policy()(wls_state,fomObj_);
  };
};
}}}





