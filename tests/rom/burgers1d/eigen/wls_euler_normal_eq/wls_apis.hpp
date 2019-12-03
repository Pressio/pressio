// this is the same as in other file
#include "wls_policies.hpp"
namespace pressio{ namespace rom{ namespace wls{
template<typename fom_type, typename wls_state_type, typename decoder_t, typename fom_state_container_t, typename hessian_gradient_pol_t>
class WlsSystemJTJApi{
private:
  const hessian_gradient_pol_t & hessian_gradient_polObj_;
//  const fom_type & fomObj_;
public:
  using scalar_t	= typename fom_type::scalar_type;
  using scalar_type = scalar_t;
  using fom_native_state_t      = typename fom_type::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;

  fom_type & appObj_;
  

  using hessian_type = typename hessian_gradient_pol_t::hess_t;
  using gradient_type = typename hessian_gradient_pol_t::gradient_t;
  using state_type = gradient_type;
  using fom_state_reconstr_t    = pressio::rom::FomStateReconstructor<fom_state_t, decoder_t>;
  using decoder_jac_t = typename decoder_t::jacobian_t;
  

  fom_state_reconstr_t  fomStateReconstructorObj_;
  fom_state_reconstr_t & fomStateReconstructorObjPtr;
  mutable fom_state_container_t  fomStateContainerObj;
  //fom_state_container_t & fomStateContainerObjPtr;
  scalar_t dt;
  int romSize; 
  int fomSize;
  int numStepsInWindow; 
  WlsSystemJTJApi(fom_type & appObj, const hessian_gradient_pol_t & hessian_gradient_polObj, const decoder_t & decoderObj, fom_state_t & yFOM, scalar_t dt, int numStepsInWindow, int romSize, int fomSize) : appObj_(appObj), hessian_gradient_polObj_(hessian_gradient_polObj),fomStateReconstructorObj_(yFOM,decoderObj),fomStateReconstructorObjPtr(fomStateReconstructorObj_),fomStateContainerObj(2,yFOM)//,fomStateContainerObjPtr(fomStateContainerObj_)
{
  this->dt = dt;
  this->numStepsInWindow = numStepsInWindow;
  this->romSize = romSize;
  this->fomSize = fomSize;
//  this->fomStateContainerObj_(2,yFOM);
 
}

  void computeHessianAndGradient(wls_state_type & wls_state, wls_state_type & wls_state_ic, hessian_type & hessian, gradient_type & gradient, const pressio::solvers::Norm & normType  = ::pressio::solvers::Norm::L2, scalar_type rnorm=0.) const{
    hessian_gradient_polObj_(appObj_,wls_state,wls_state_ic,hessian,gradient,fomStateContainerObj,fomStateReconstructorObjPtr,dt,numStepsInWindow);
  }

  hessian_type createHessianObject(const wls_state_type & stateIn) const{
    hessian_type H(romSize*numStepsInWindow,romSize*numStepsInWindow);  //how do we do this for arbitrary types?
    return H;} 

  gradient_type createGradientObject(const wls_state_type & stateIn) const{
    gradient_type g(romSize*numStepsInWindow);  //how do we do this for arbitrary types?
    return g;} 


};





/*

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
*/
}}}





