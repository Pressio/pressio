// this is the same as in other file
#include "wls_policies.hpp"
namespace pressio{ namespace rom{ namespace wls{
template<typename fom_type, 
	typename wls_state_type, 
	typename decoder_t, 
	typename hessian_gradient_pol_t, 
	typename aux_states_container_t, 
	typename hessian_t>
class WlsSystemHessianAndGradientApi{
public:
  using gradient_type = wls_state_type;
  using scalar_t	= typename fom_type::scalar_type;
  using scalar_type = scalar_t;
  using fom_native_state_t      = typename fom_type::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using state_type = wls_state_type;
  using hessian_type = hessian_t;
  using fom_state_reconstr_t    = pressio::rom::FomStateReconstructor<fom_state_t, decoder_t>;
  using decoder_jac_t = typename decoder_t::jacobian_t;
private:
  const hessian_gradient_pol_t & hessian_gradient_polObj_;
  fom_type & appObj_;
  fom_state_reconstr_t  fomStateReconstructorObj_;
  mutable aux_states_container_t  auxStatesContainer;
  mutable wls_state_type wlsStateIC;
  scalar_t dt;
  mutable scalar_t ts = 0.;
  int romSize; 
  int fomSize;
  int numStepsInWindow;
  int time_stencil_size;
  int windowNum = 0;
public: 
  WlsSystemHessianAndGradientApi(fom_type & appObj, 
				const hessian_gradient_pol_t & hessian_gradient_polObj, 
				const decoder_t & decoderObj,
				fom_state_t & yFOM, 
				scalar_t dt, 
				int numStepsInWindow, 
				int romSize, 
				int fomSize, 
				const int time_stencil_size) : 
				appObj_(appObj), 
				hessian_gradient_polObj_(hessian_gradient_polObj),
				fomStateReconstructorObj_(yFOM,decoderObj),
				wlsStateIC(romSize*(time_stencil_size-1)), 
				auxStatesContainer(yFOM){
    this->dt = dt;
    this->numStepsInWindow = numStepsInWindow;
    this->romSize = romSize;
    this->fomSize = fomSize;
    this->time_stencil_size = time_stencil_size;
    this->wlsStateIC.setZero();
  }
  void computeHessianAndGradient(const state_type & wls_state,hessian_type & hessian, gradient_type & gradient, const pressio::solvers::Norm & normType  = ::pressio::solvers::Norm::L2, scalar_type rnorm=0.) const{
    hessian_gradient_polObj_(appObj_,wls_state,wlsStateIC,hessian,gradient,fomStateReconstructorObj_,dt,numStepsInWindow,this->ts,auxStatesContainer);
  }
  hessian_type createHessianObject(const wls_state_type & stateIn) const{
    hessian_type H(romSize*numStepsInWindow,romSize*numStepsInWindow);  //how do we do this for arbitrary types?
    return H;} 
  gradient_type createGradientObject(const wls_state_type & stateIn) const{
    gradient_type g(romSize*numStepsInWindow);  //how do we do this for arbitrary types?
    return g;} 



  template <typename solverType>
  void advanceOneWindow(wls_state_type & wlsState, solverType & solver, int  windowNum) {
    this->windowNum = windowNum;
    this->ts  = windowNum*this->dt*this->numStepsInWindow; 
    //solver.my_gauss_newton(*this,wlsState,romSize,numStepsInWindow);
    solver.solve(*this,wlsState);
    //sbar = std::min(this->time_stencil_size-1,numStepsInWindow);
    int start = std::max(0,time_stencil_size - 1 - numStepsInWindow);
    for (int i = 0; i < start; i++){
      auto wlsTmpState = ::pressio::containers::span(wlsStateIC,i*romSize,romSize);
      auto wlsTmpState2 = ::pressio::containers::span(wlsStateIC,(i+1)*romSize,romSize);
      for (int k=0; k < romSize; k++){
        wlsTmpState[k] = wlsTmpState2[k];}}
    for (int i = start ; i < this->time_stencil_size-1; i++){
      auto wlsTmpState = ::pressio::containers::span(wlsStateIC,i*romSize,romSize);
      auto wlsTmpState2 = ::pressio::containers::span(wlsState,( numStepsInWindow - this->time_stencil_size + 1 + i )*romSize,romSize);
      for (int k=0; k < romSize; k++){
        wlsTmpState[k] = wlsTmpState2[k];}}
      //(*(this->wlsStateIC).data()).block(0,i,romSize,1) = (*wlsState.data()).block(0,numStepsInWindow - this->time_stencil_size + 1 + i , romSize ,1);}
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
}}}





