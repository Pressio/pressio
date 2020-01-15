#include "../../ode/select_ode_helper.hpp"
#include "../../policies/hess_and_grad_policies/hess_and_grad_policies.hpp"

#include "SOLVERS_LINEAR"

namespace pressio{ namespace rom{ namespace wls{ namespace impl{
template<typename fom_type,        
	typename wls_state_type, 
	typename decoder_t,
  typename ode_tag, 
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
  using hessian_gradient_pol_t = ::pressio::rom::wls::hessian_gradient_policy<fom_type,decoder_t>;

  // Assertions 
  static_assert(std::is_same<wls_state_type, ::pressio::containers::Vector<fom_native_state_t> >::value, "WLS_state_type  must be of type ::pressio::containers::Vector");
  static_assert(std::is_same<hessian_type,pressio::containers::Matrix<Eigen::Matrix<scalar_t, -1, -1>> >::value, "hessian_type must be of type pressio::containers::Matrix<Eigen::Matrix<scalar_t, -1, -1>>");

private:
  const hessian_gradient_pol_t hessian_gradient_polObj_; 
  const fom_type & appObj_;
  const fom_state_reconstr_t  fomStateReconstructorObj_;

  //create an ode object that contains information on stencil width, time discrete resiudal, time discrete jacobian, etc. 
  using ode_policies_t = ::pressio::rom::wls::ode::helpers::ode_policies_t<ode_tag,fom_state_t, wls_state_type>; 

  scalar_t dt_;
  scalar_t ts_ = 0.;      //time at start of window
  int romSize_;           //size of generalized coordinates 
  int numStepsInWindow_;  //number of discrete time instances in a window
  int timeStencilSize_;   //will serve as an alias to the ode policiy
  int windowNum_ = 0;     // active window index
  int step_s_ = 0;        // global step number
public: 
  ode_policies_t odeObj_;
  WlsSystemHessianAndGradientApi(const fom_type & appObj, 
                                 const decoder_t & decoderObj,
                                 const fom_state_t & yFOM_IC, 
                                 const fom_state_t & yFOM_Ref, 
                                 const int numStepsInWindow)
                                 : 
                                 romSize_(decoderObj.getReferenceToJacobian().numVectors()),
                                 odeObj_(romSize_,yFOM_IC),
                                 appObj_(appObj),
                                 fomStateReconstructorObj_(yFOM_Ref,decoderObj),
                                 hessian_gradient_polObj_( appObj, yFOM_IC,numStepsInWindow,odeObj_.state_stencil_size_,decoderObj)
                                 {
                                   this->numStepsInWindow_ = numStepsInWindow; //set for convenience
                                   this->timeStencilSize_ = odeObj_.state_stencil_size_; //set for convenience
                                   // Set initial condition based on L^2 projection onto trial space
                                   // note that wlsStateIC_[-romSize:end] contains nm1, wlsStateIC[-2*romSize:-romSize] contains nm2 entry, etc.
                                   auto wlsInitialStateNm1 = ::pressio::containers::span(odeObj_.wlsStateIC_,romSize_*(odeObj_.state_stencil_size_-2),romSize_);
                                   initializeCoeffs(decoderObj.getReferenceToJacobian(), wlsInitialStateNm1 , yFOM_IC,yFOM_Ref);
                                 }


  // Function to initialize the coefficients for a given FOM IC and reference state
  // computes the ICs via optimal L^2 projection, phi^T phi xhat = phi^T(x - xRef)
  template <typename basis_t, typename wls_stateview_t>
  void initializeCoeffs(const basis_t & phi, wls_stateview_t & wlsStateIC, const fom_state_t & yFOM_IC,const fom_state_t & yRef ){
    using solver_tag   = pressio::solvers::linear::direct::ColPivHouseholderQR;
    using linear_solver_t = pressio::solvers::direct::EigenDirect<solver_tag, hessian_t>;
    linear_solver_t linear_solver; //initialize linear solver. This is only used here 
    hessian_t H(this->romSize_,this->romSize_); // create the system matrix, phi^T phi
    ::pressio::containers::ops::dot(phi,phi,H);
    fom_state_t b(yFOM_IC); //create a vector to store yFOM - yRef.
    pressio::containers::ops::do_update(b,1.,yRef,-1.);
    auto r = pressio::containers::ops::dot(phi,b); // compute phi^T b
    wls_state_type y_l2(this->romSize_); //currently have no way of passing a span to the linear solver
    linear_solver.solveAllowMatOverwrite(H, r, y_l2); //solve system for optimal L2 projection
    for (int i=0; i< this->romSize_; i++){
      wlsStateIC[i] = y_l2[i]; //set ICs
    }
  }

  // Policity to compute Hessian and gradient
  void computeHessianAndGradient(const state_type & wls_state,
                                 hessian_type & hessian, 
                                 gradient_type & gradient, 
                                 const pressio::solvers::Norm & normType  = ::pressio::solvers::Norm::L2, 
                                 scalar_type & rnorm=0.) const{
    rnorm = 0.;
    hessian_gradient_polObj_(this->appObj_,
                             this->odeObj_,
                             wls_state,
                             this->odeObj_.wlsStateIC_,
                             hessian,
                             gradient,
                             this->fomStateReconstructorObj_,
                             this->dt_,
                             this->numStepsInWindow_,
                             this->ts_,
                             this->step_s_,
                             rnorm);
  }

  hessian_type createHessianObject(const wls_state_type & stateIn) const{
    hessian_type H(romSize_*numStepsInWindow_,romSize_*numStepsInWindow_);  //how do we do this for arbitrary types?
    return H;} 
  gradient_type createGradientObject(const wls_state_type & stateIn) const{
    gradient_type g(romSize_*numStepsInWindow_);  //how do we do this for arbitrary types?
    return g;} 



  // method to advance one window. We may want to put this into some type of window stepper class 
  // if we want to have more complex stepping  
  template <typename solverType>
  void advanceOneWindow(wls_state_type & wlsState,solverType & solver,const int  & windowNum, scalar_t dt) {
    this->dt_ = dt; //set time step
    this->windowNum_ = windowNum; //set window number
    this->ts_  = windowNum*this->dt_*this->numStepsInWindow_;  //set starting time
    this->step_s_  = windowNum*this->numStepsInWindow_;  //set step number
    solver.solve(*this,wlsState); //solve system
    // Loop to update the the wlsStateIC vector. 
    // If we add multistep explicit methods, need to add things here.
    int start = std::max(0,this->timeStencilSize_ - 1 - this->numStepsInWindow_); 
    for (int i = 0; i < start; i++){
      auto wlsTmpState = ::pressio::containers::span(this->odeObj_.wlsStateIC_,i*this->romSize_,this->romSize_);
      auto wlsTmpState2 = ::pressio::containers::span(this->odeObj_.wlsStateIC_,(i+1)*this->romSize_,this->romSize_); 
      for (int k=0; k < romSize_; k++){
        wlsTmpState[k] = wlsTmpState2[k];}}
    for (int i = start ; i < this->timeStencilSize_-1; i++){
      auto wlsTmpState = ::pressio::containers::span(this->odeObj_.wlsStateIC_,i*this->romSize_,this->romSize_);
      auto wlsTmpState2 = ::pressio::containers::span(wlsState,( numStepsInWindow_ - this->timeStencilSize_ + 1 + i )*this->romSize_,this->romSize_);
      for (int k=0; k < romSize_; k++){
        wlsTmpState[k] = wlsTmpState2[k];}}
    std::cout << " Window " << windowNum << " completed " << std::endl;
  }

};

}}}}





