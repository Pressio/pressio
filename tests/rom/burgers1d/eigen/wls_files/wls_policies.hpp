#include "wls_local_policies.hpp"
#include "wls_functions.hpp"


/* 
This header file contains policies for computing the hessian and gradients in the WLS system.
The policies are templated on the policies for compputing the local residuals and jacobians as defined in the wls_local_policies.hpp file
*/

namespace pressio{ namespace rom{ namespace wls{ namespace impl{
template< typename wls_state_type, typename fom_type, typename hess_type, typename gradient_type, typename decoder_t,typename local_residual_policy_t,typename local_jacobian_policy_t>
struct hessian_gradient_policy{
public:
  using hess_t = hess_type;
  using gradient_t = gradient_type;
  using fom_native_state_t      = typename fom_type::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using fom_state_container_t = std::vector<fom_state_t>;
  using fom_state_reconstr_t    = pressio::rom::FomStateReconstructor<fom_state_t, decoder_t>;
  using scalar_t	= typename fom_type::scalar_type;
  using basis_t  = typename decoder_t::jacobian_t;
  using residual_t = fom_state_t; //maybe we need to change this for hyper reduction?
  // Create types for Jacobians and Jacobian containers
  using jac_t = typename pressio::containers::Matrix<Eigen::Matrix<scalar_t, -1, -1>>;
  using wls_jacs_t = std::vector<jac_t>;

  // Initialize required objects
  int romSize;
  int fomSize;
  int time_stencil_size;
  basis_t & phi_;
  mutable residual_t residual;
  mutable wls_jacs_t wlsJacs;
  local_residual_policy_t residualPolicy;
  local_jacobian_policy_t jacobianPolicy;

  hessian_gradient_policy(int romSize,  int fomSize, int numStepsInWindow, int time_stencil_size,  basis_t & phi) : wlsJacs(time_stencil_size,jac_t(fomSize,romSize)), phi_(phi), residual(fomSize){ 
    this->romSize = romSize;    
    this->time_stencil_size = time_stencil_size;
    this->fomSize = fomSize;
  }

  void operator()(fom_type & appObj, const wls_state_type  & wlsState, wls_state_type & wlsStateIC, hess_type & hess, gradient_type & gradient,fom_state_container_t &  fomStateContainerObj_, fom_state_reconstr_t & fomStateReconstrObj_,scalar_t dt, std::size_t numStepsInWindow) const{
  int n = 0;
  const auto wlsCurrentState = wlsState.viewColumnVector(n);
  fomStateReconstrObj_(wlsCurrentState,fomStateContainerObj_[1]);
  const auto wlsInitialState = wlsStateIC.viewColumnVector(0);
  fomStateReconstrObj_(wlsInitialState,fomStateContainerObj_[0]);
  residualPolicy(appObj,residual,fomStateContainerObj_,n*dt,dt);
  for (int i = 0; i < time_stencil_size; i++){
    jacobianPolicy(appObj,wlsJacs[time_stencil_size - i - 1],phi_,fomStateContainerObj_,n*dt,dt,i);
  }
  local_mat_update(wlsJacs[1],wlsJacs[1],hess,n*romSize,n*romSize,romSize,romSize);
  //(*hess.data()).block(n*romSize,n*romSize,romSize,romSize) = (*wlsJacs[1].data()).transpose() * (*wlsJacs[1].data());
  (*gradient.data()).block(n*romSize,0,romSize,1) = -(*wlsJacs[1].data()).transpose()*(*residual.data());
  for (int n = 1; n < numStepsInWindow; n++){
    // === reconstruct FOM states ========
    *fomStateContainerObj_[0].data() = *fomStateContainerObj_[1].data();
    const auto wlsCurrentState = wlsState.viewColumnVector(n);
    fomStateReconstrObj_(wlsCurrentState,fomStateContainerObj_[1]);
    // == Evaluate residual ============
    residualPolicy(appObj,residual,fomStateContainerObj_,n*dt,dt);
    // == Evaluate Jacobian (Currently have jacobian w.r.p to previous state hard coded outside of loop)
    jacobianPolicy(appObj,wlsJacs[1],phi_,fomStateContainerObj_,n*dt,dt,0);
    (*gradient.data()).block(n*romSize,0,romSize,1) = -(*wlsJacs[1].data()).transpose()*(*residual.data());
    for (int i=1; i < time_stencil_size; i++){
      (*gradient.data()).block((n-i)*romSize,0,romSize,1) -= (*wlsJacs[time_stencil_size-i-1].data()).transpose()*(*residual.data());}
    // == Assembel local component of global Hessian //
    for (int i=0; i < time_stencil_size; i++){
      for (int j=0; j <= i; j++){
         local_mat_update(wlsJacs[time_stencil_size-i-1],wlsJacs[time_stencil_size-j-1],hess,(n-i)*romSize,(n-j)*romSize,romSize,romSize);
         //(*hess.data()).block((n-i)*romSize,(n-j)*romSize,romSize,romSize) += (*wlsJacs[time_stencil_size-i-1].data()).transpose() * (*wlsJacs[time_stencil_size-j-1].data());
         block_sym_update(hess,(n-i)*romSize,(n-j)*romSize,romSize,romSize); 
//         if (i != j){
//          (*hess.data()).block((n-j)*romSize,(n-i)*romSize,romSize,romSize) =  (*hess.data()).block((n-i)*romSize,(n-j)*romSize,romSize,romSize).transpose();}
        }
    }
  } 
}
};

/*
template< typename wls_state_type, typename fom_type, typename residual_type, typename decoder_t>
struct residual_policy_naive{
public:
  using residual_t = residual_type;
  using fom_native_state_t      = typename fom_type::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using fom_state_reconstr_t    = pressio::rom::FomStateReconstructor<fom_state_t, decoder_t>;
  using fom_state_container_t = typename pressio::ode::AuxStatesContainer<fom_state_t, 2>;
  using scalar_t	= typename fom_type::scalar_type;

  void operator()(fom_type & appObj, wls_state_type  & wlsState, residual_t & residual, fom_state_container_t &  fomStateContainerObj_, fom_state_reconstr_t & fomStateReconstrObj_,scalar_t dt) const { 
  std::cout << "Computed the (void) residual through the naive policy" << std::endl;
  for (int i=1; i< 5; i++){
    fomStateReconstrObj_(wlsState[i],fomStateContainerObj_[0]);
    fomStateReconstrObj_(wlsState[i-1],fomStateContainerObj_[1]);
    std::cout << fomStateContainerObj_[0][i]<< std::endl;
    appObj.timeDiscreteResidual(i,dt*i,dt,*fomStateContainerObj_[0].data(),*fomStateContainerObj_[0].data(),*fomStateContainerObj_[1].data());
    }
  }
 
  residual_t evaluate(wls_state_type  & wls_state) const { 
  std::cout << "Computed the residual through the naive policy" << std::endl; 
  residual_t a;// placeholder to return residual
  return a;
}
};


template< typename wls_state_type, typename fom_type, typename jacobian_type>
struct jacobian_policy_naive{
public:
  using jacobian_t = jacobian_type;
  void operator()(wls_state_type  & wlsState, jacobian_t & jacobian) const { 
  std::cout << "Computed the (void) jacobian through the naive policy" << std::endl; 
}
};

*/


}}}}

