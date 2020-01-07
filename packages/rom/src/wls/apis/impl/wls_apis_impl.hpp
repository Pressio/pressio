// this is the same as in other file
#include "../../policies/wls_policies.hpp"

namespace pressio{ namespace rom{ namespace wls{ namespace impl{
template<typename fom_type, 
	typename wls_state_type, 
	typename decoder_t,
        typename ode_t, 
	typename hessian_gradient_pol_t, 
	typename aux_states_container_t, 
	typename hessian_t>
class WlsSystemHessianAndGradientApi{
public:
  using gradient_type = wls_state_type;
  using scalar_t	= typename fom_type::scalar_type;
  using scalar_type = scalar_t;//
  using fom_native_state_t      = typename fom_type::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using state_type = wls_state_type;
  using hessian_type = hessian_t;
  using fom_state_reconstr_t    = pressio::rom::FomStateReconstructor<fom_state_t, decoder_t>;
  using decoder_jac_t = typename decoder_t::jacobian_t;
private:
  const hessian_gradient_pol_t hessian_gradient_polObj_; 
  const fom_type & appObj_;
  const fom_state_reconstr_t  fomStateReconstructorObj_;
  mutable aux_states_container_t  auxStatesContainer_; //understand mutable keyword better
  wls_state_type wlsStateIC_;
  scalar_t dt_;
  scalar_t ts_ = 0.;
  int romSize_; 
  int fomSize_;
  int numStepsInWindow_;
  int timeStencilSize_;
  int windowNum_ = 0;
  int step_s_ = 0;
public: 
  WlsSystemHessianAndGradientApi(const fom_type & appObj, 
                                 const decoder_t & decoderObj,
                                 fom_state_t & yFOM, 
                                 const int numStepsInWindow, 
                                 const int timeStencilSize) 
				 : 
                                 appObj_(appObj),
                                 fomSize_(yFOM.size()),
                                 romSize_((decoderObj.getReferenceToJacobian()).numVectors()),
                                 fomStateReconstructorObj_(yFOM,decoderObj),
                                 wlsStateIC_( ( (decoderObj.getReferenceToJacobian()).numVectors() ) *(timeStencilSize-1)),
                                 auxStatesContainer_(yFOM),
                                 hessian_gradient_polObj_( ( (decoderObj.getReferenceToJacobian()).numVectors() ) ,yFOM.size(),numStepsInWindow,timeStencilSize,decoderObj)
				 {
                                   this->numStepsInWindow_ = numStepsInWindow;
                                   this->timeStencilSize_ = timeStencilSize;
                                   this->wlsStateIC_.setZero();
  			         }
  void computeHessianAndGradient(const state_type & wls_state,
                                 hessian_type & hessian, 
                                 gradient_type & gradient, 
                                 const pressio::solvers::Norm & normType  = ::pressio::solvers::Norm::L2, 
                                 scalar_type rnorm=0.) const{
    hessian_gradient_polObj_(this->appObj_,
                             wls_state,
                             this->wlsStateIC_,
                             hessian,
                             gradient,
                             this->fomStateReconstructorObj_,
                             this->dt_,
                             this->numStepsInWindow_,
                             this->ts_,
                             this->auxStatesContainer_,
                             this->step_s_);
  }
  //We have to put static assertion for gradient and hessian to be pressio wrappers for specific types
  hessian_type createHessianObject(const wls_state_type & stateIn) const{
    hessian_type H(romSize_*numStepsInWindow_,romSize_*numStepsInWindow_);  //how do we do this for arbitrary types?
    return H;} 
  gradient_type createGradientObject(const wls_state_type & stateIn) const{
    gradient_type g(romSize_*numStepsInWindow_);  //how do we do this for arbitrary types?
    return g;} 



  template <typename solverType>
  void advanceOneWindow(wls_state_type & wlsState,solverType & solver,const int  & windowNum, scalar_t dt) {
    this->dt_ = dt;
    this->windowNum_ = windowNum;
    this->ts_  = windowNum*this->dt_*this->numStepsInWindow_; 
    this->step_s_  = windowNum*this->numStepsInWindow_; 
    solver.solve(*this,wlsState);
    int start = std::max(0,timeStencilSize_ - 1 - numStepsInWindow_); 
    for (int i = 0; i < start; i++){
      auto wlsTmpState = ::pressio::containers::span(wlsStateIC_,i*romSize_,romSize_);
      auto wlsTmpState2 = ::pressio::containers::span(wlsStateIC_,(i+1)*romSize_,romSize_); 
      for (int k=0; k < romSize_; k++){
        wlsTmpState[k] = wlsTmpState2[k];}}
    for (int i = start ; i < this->timeStencilSize_-1; i++){
      auto wlsTmpState = ::pressio::containers::span(wlsStateIC_,i*romSize_,romSize_);
      auto wlsTmpState2 = ::pressio::containers::span(wlsState,( numStepsInWindow_ - this->timeStencilSize_ + 1 + i )*romSize_,romSize_);
      for (int k=0; k < romSize_; k++){
        wlsTmpState[k] = wlsTmpState2[k];}}
    std::cout << " Window " << windowNum << " completed " << std::endl;
  }

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
}}}}





