#include "../../../ode/select_ode_helper.hpp"
/* 
This header file contains the class object used for computing the hessian and gradients in the WLS system.
This class creates an ode object from ::pressio::rom::wls::ode , which then contains routines for the residual, jacobian, etc.

The hessian_gradient policy is responsible for assembling Jw^T Jw and J^T rw, where Jw is the reduced windowed Jacobian and rw the residual. 
The policy initializes in memory:
wlsJacs: this is a vector container of n_s -1 local Jacobians, J, where n_s is the width of the time stencil (e.g., n_s = 3 for BDF2)
  residual: this is the residual vector for the FOM 
  yFOM_current: this is a working variable for the fom state.
*/
namespace pressio{ namespace rom{ namespace wls{ namespace impl{
template<typename fom_type,
	 typename decoder_t,
	 typename ode_tag>
class hessian_gradient_policy{
private:
  using fom_native_state_t      = typename fom_type::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using scalar_t                = typename fom_type::scalar_type;
  using jacobian_t              = typename decoder_t::jacobian_t; 
  using jac_t                   =  jacobian_t;
  using residual_t              = fom_state_t; //maybe we need to change this for hyper reduction?
  // Create types for Jacobians and Jacobian containers
  using wls_jacs_t = std::vector<jac_t>; //
  using nm1 = ::pressio::ode::nMinusOne;
  using nm2 = ::pressio::ode::nMinusTwo;
  using ode_policies_t = ::pressio::rom::wls::ode::helpers::ode_policies_t<ode_tag,fom_state_t>;
  ode_policies_t odeObj_;
  mutable wls_jacs_t wlsJacs_;
  mutable fom_state_t yFOM_current_;
  mutable residual_t residual_;
  int romSize_;
  int fomSize_;
  int time_stencil_size;
  const jac_t & phi_; 
  //==========================================
  // Function to do c += A^T b
  template <typename Mat, 
            typename Vec>
  void local_vec_update(const Mat & A,
                        const  Vec & b, 
                        Vec  & cOut,
                        const int & scol,
                        const int & colSize) const{
    auto gradientView = ::pressio::containers::span(cOut,scol,colSize);
    auto tmp = pressio::containers::ops::dot(A, b);
    for (int k =0; k< colSize; k++){
      gradientView[k] += tmp[k];}
    }
public:

  hessian_gradient_policy(const fom_type & appObj,
                          const fom_state_t & yFOM,
                          int numStepsInWindow, 
                          int time_stencil_size,  
                          const decoder_t & decoderObj) : 
                          wlsJacs_(time_stencil_size,decoderObj.getReferenceToJacobian()),  // construct Jacobians off of decoder type. // jac_t(fomSize,romSize)),
                          phi_(decoderObj.getReferenceToJacobian()), 
                          residual_( appObj.velocity( *yFOM.data() , ::pressio::utils::constants::zero<scalar_t>()) ),
                          yFOM_current_(yFOM),
                          odeObj_(romSize_) 
    { 
    this->romSize_ = decoderObj.getReferenceToJacobian().numVectors(); 
    this->time_stencil_size = time_stencil_size;
    this->fomSize_ = yFOM.size();
    }
  template <typename wls_state_type,
            typename fom_state_reconstr_t, 
            typename aux_states_container_t, 
            typename hess_type, 
            typename gradient_type>
  void operator()(const fom_type & appObj, 
                  const wls_state_type  & wlsState,
                  const wls_state_type & wlsStateIC, 
                  hess_type & hess, 
                  gradient_type & gradient, 
                  const fom_state_reconstr_t & fomStateReconstrObj_,
                  const scalar_t dt, 
                  const std::size_t numStepsInWindow, 
                  const scalar_t ts , 
                  aux_states_container_t & auxStatesContainer, 
                  const int step_s,
                  scalar_t & rnorm) const{
  int n = 0;
  scalar_t t = ts + n*dt;
  int step = step_s + n;
  hess.setZero(); //ensure zero as we use += 
  gradient.setZero();

  const auto wlsCurrentState = ::pressio::containers::span(wlsState,0,romSize_); //get access to the state at the first window
  fomStateReconstrObj_(wlsCurrentState,yFOM_current_);
  odeObj_.updateStatesFirstStep(wlsStateIC,fomStateReconstrObj_,auxStatesContainer); // reconstruct the FOM states from the previous window/ICs 

  odeObj_.time_discrete_residual(appObj,yFOM_current_,residual_,auxStatesContainer,ts,dt,step); //compute the time discrete residual
  rnorm += ::pressio::containers::ops::norm2(residual_); //increment the norm

  hess_type C(romSize_,romSize_); // this is a temporary work around for the += issue.
  for (int i = 0; i < time_stencil_size; i++){
    odeObj_.time_discrete_jacobian(appObj,yFOM_current_,wlsJacs_[time_stencil_size - i - 1],phi_,auxStatesContainer,n*dt,dt,step,i); //compute the Jacobian for each state in the stencil
  }
  // add to local block of hessian
  auto hess_block = ::pressio::containers::subspan( hess,std::make_pair( n*romSize_,(n+1)*romSize_ ) , std::make_pair( n*romSize_,(n+1)*romSize_) );
  ::pressio::containers::ops::dot(wlsJacs_[time_stencil_size-1],wlsJacs_[time_stencil_size-1], hess_block);

  local_vec_update(wlsJacs_[time_stencil_size-1],residual_,gradient,n*romSize_,romSize_);  // compute gradient[n*romSize_:(n+1)*romSize] += J^T r
  for (int n = 1; n < numStepsInWindow; n++){
    // === reconstruct FOM states ========
    odeObj_.updateStatesNStep(yFOM_current_,auxStatesContainer); 

    const auto wlsCurrentState = ::pressio::containers::span(wlsState,n*romSize_,romSize_);
    fomStateReconstrObj_(wlsCurrentState,yFOM_current_);
    // == Evaluate residual ============
    t = ts + n*dt;
    step = step_s + n; 
    odeObj_.time_discrete_residual(appObj,yFOM_current_,residual_,auxStatesContainer,t,dt,step);
    rnorm += ::pressio::containers::ops::norm2(residual_);
    for (int i = 0; i < time_stencil_size; i++){
      odeObj_.time_discrete_jacobian(appObj,yFOM_current_,wlsJacs_[time_stencil_size - i - 1],phi_,auxStatesContainer,n*dt,dt,step,i); //compute the Jacobian for each state in the stencil
    }
    // == Update everything
    int sbar = std::min(n,time_stencil_size);
    for (int i=0; i < sbar; i++){
      local_vec_update(wlsJacs_[time_stencil_size-i-1],residual_,gradient,(n-i)*romSize_,romSize_);}
    // == Assemble local component of global Hessian //
    for (int i=0; i < sbar; i++){
      for (int j=0; j <= i; j++){
         auto hess_block = ::pressio::containers::subspan(hess, std::make_pair( (n-i)*romSize_,(n-i+1)*romSize_ ) , std::make_pair( (n-j)*romSize_,(n-j+1)*romSize_ ) );
         ::pressio::containers::ops::dot(wlsJacs_[time_stencil_size-i-1],wlsJacs_[time_stencil_size-j-1], C);
         for (int k = 0; k<romSize_; k++){
           for (int l = 0; l<romSize_; l++){
             hess_block(k,l) = hess_block(k,l) + C(k,l); 
            }
         }
         auto hess_block2 = ::pressio::containers::subspan(hess, std::make_pair( (n-j)*romSize_,(n-j+1)*romSize_) , std::make_pair( (n-i)*romSize_,(n-i+1)*romSize_) );
         for (int k = 0; k<romSize_; k++){
           for (int l = 0; l<romSize_; l++){
             hess_block2(l,k) = hess_block(k,l); 
            }
         }
      }
    }
  } 
}


  template <typename wls_state_type,
            typename fom_state_reconstr_t, 
            typename aux_states_container_t> 
  void computeResidualNorm(const fom_type & appObj, 
                  const wls_state_type  & wlsState,
                  const wls_state_type & wlsStateIC, 
                  const fom_state_reconstr_t & fomStateReconstrObj_,
                  const scalar_t dt, 
                  const std::size_t numStepsInWindow, 
                  const scalar_t ts , 
                  aux_states_container_t & auxStatesContainer, 
                  const int step_s,
                  scalar_t & rnorm) const{
  rnorm = ::pressio::utils::constants::zero<scalar_t>();
  int n = 0;
  scalar_t t = ts + n*dt;
  int step = step_s + n;
  const auto wlsCurrentState = ::pressio::containers::span(wlsState,0,romSize_);
  fomStateReconstrObj_(wlsCurrentState,yFOM_current_);
  odeObj_.updateStatesFirstStep(wlsStateIC,fomStateReconstrObj_,auxStatesContainer);
  odeObj_.time_discrete_residual(appObj,yFOM_current_,residual_,auxStatesContainer,ts,dt,step);
  rnorm += ::pressio::containers::ops::norm2(residual_);
  for (int n = 1; n < numStepsInWindow; n++){
    // === reconstruct FOM states ========
    odeObj_.updateStatesNStep(yFOM_current_,auxStatesContainer); 
    const auto wlsCurrentState = ::pressio::containers::span(wlsState,n*romSize_,romSize_);
    fomStateReconstrObj_(wlsCurrentState,yFOM_current_);
    // == Evaluate residual ============
    t = ts + n*dt;
    step = step_s + n; 
    odeObj_.time_discrete_residual(appObj,yFOM_current_,residual_,auxStatesContainer,t,dt,step);
    rnorm += ::pressio::containers::ops::norm2(residual_);
  }
}


};



}}}}

