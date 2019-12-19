#include "wls_local_policies.hpp"
//#include "wls_functions.hpp"


/* 
This header file contains policies for computing the hessian and gradients in the WLS system.
The policies are templated on the policies for compputing the local residuals and jacobians as defined in the wls_local_policies.hpp file
*/

/*
The hessian_gradient policy is responsible for assembling Jw^T Jw and J^T rw, where Jw is the reduced windowed Jacobian and rw the residual. 
The policy initializes in memory:
	wlsJacs: this is a vector container of n_s -1 local Jacobians, J, where n_s is the width of the time stencil (e.g., n_s = 3 for BDF2)
  residual: this is the residual vector for the FOM 
  yFOM_current: this is a working variable for the fom state.
*/
namespace pressio{ namespace rom{ namespace wls{ namespace impl{
template<typename fom_type, typename jacobian_t,typename local_residual_policy_t,typename local_jacobian_policy_t, typename  ode_tag>
class hessian_gradient_policy{
private:
  using fom_native_state_t      = typename fom_type::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  //using fom_state_reconstr_t    = pressio::rom::FomStateReconstructor<fom_state_t, decoder_t>;
  using scalar_t	= typename fom_type::scalar_type;
  using jac_t  =  jacobian_t;
  using residual_t = fom_state_t; //maybe we need to change this for hyper reduction?
  // Create types for Jacobians and Jacobian containers
  using wls_jacs_t = std::vector<jac_t>;
  using nm1 = ode::nMinusOne;
  using nm2 = ode::nMinusTwo;
  mutable wls_jacs_t wlsJacs;
  mutable fom_state_t yFOM_current;
  int romSize;
  int fomSize;
  int time_stencil_size;
  jac_t & phi_;
  mutable residual_t residual;
  local_residual_policy_t residualPolicy;
  local_jacobian_policy_t jacobianPolicy;
  ode_tag ode;
  //=============================================
  // Functions to update the state in the WLS loop
  // BDF 2
  template <typename wls_state_type,typename fom_state_reconstr_t,  typename aux_states_container_t> 
  void updateStatesFirstStep(const wls_state_type & wlsStateIC, const fom_state_reconstr_t & fomStateReconstr, aux_states_container_t & auxStatesContainer, ::pressio::ode::implicitmethods::BDF2  tag) const{
    const auto wlsInitialStateNm2 = ::pressio::containers::span(wlsStateIC,0,romSize);
    const auto wlsInitialStateNm1 = ::pressio::containers::span(wlsStateIC,romSize,romSize);
    auto & fomStateNm2 = auxStatesContainer.get(nm2());
    auto & fomStateNm1 = auxStatesContainer.get(nm1());
    fomStateReconstr(wlsInitialStateNm2,fomStateNm2);
    fomStateReconstr(wlsInitialStateNm1,fomStateNm1);
  }

  template <typename aux_states_container_t> 
  void updateStatesNStep(const fom_state_t & yFOM_current, aux_states_container_t & auxStatesContainer, const ::pressio::ode::implicitmethods::BDF2 & tag) const{
      auto & odeState_nm1 = auxStatesContainer.get(nm1());
      auto & odeState_nm2 = auxStatesContainer.get(nm2());
      ::pressio::containers::ops::deep_copy(odeState_nm1, odeState_nm2);
      ::pressio::containers::ops::deep_copy(yFOM_current, odeState_nm1);
  }
  // Euler
  template <typename wls_state_type,typename fom_state_reconstr_t, typename aux_states_container_t> 
  void updateStatesFirstStep(const wls_state_type & wlsStateIC, const fom_state_reconstr_t & fomStateReconstr, aux_states_container_t & auxStatesContainer, const ::pressio::ode::implicitmethods::Euler & tag) const{
    const auto wlsInitialStateNm1 = ::pressio::containers::span(wlsStateIC,0,romSize);
    auto & fomStateNm1 = auxStatesContainer.get(nm1());
    fomStateReconstr(wlsInitialStateNm1,fomStateNm1);
  }

  template <typename aux_states_container_t> 
  void updateStatesNStep(const fom_state_t & yFOM_current, aux_states_container_t & auxStatesContainer, const ::pressio::ode::implicitmethods::Euler  & tag) const{
      auto & odeState_nm1 = auxStatesContainer.get(nm1());
      ::pressio::containers::ops::deep_copy(yFOM_current, odeState_nm1);
  }
  //==========================================
  // Function to do c += A^T b
  template <typename Mat, typename Vec>
  void local_vec_update(const Mat & A,const  Vec & b, Vec  & cOut,const int & scol,const int & colSize) const{
    auto gradientView = ::pressio::containers::span(cOut,scol,colSize);
    auto tmp = pressio::containers::ops::dot(A, b);
    for (int k =0; k< colSize; k++){
      gradientView[k] += tmp[k];}
  }

public:
  hessian_gradient_policy(int romSize,  int fomSize, int numStepsInWindow, int time_stencil_size,  jac_t & phi) : wlsJacs(time_stencil_size,jac_t(fomSize,romSize)), phi_(phi), residual(fomSize), yFOM_current(fomSize){ 
    this->romSize = romSize;    
    this->time_stencil_size = time_stencil_size;
    this->fomSize = fomSize;
  }

  template <typename wls_state_type,typename fom_state_reconstr_t, typename  aux_states_container_t, typename hess_type, typename gradient_type>
  void operator()(const fom_type & appObj, const wls_state_type  & wlsState,const wls_state_type & wlsStateIC, hess_type & hess, gradient_type & gradient, const fom_state_reconstr_t & fomStateReconstrObj_,const scalar_t dt, const std::size_t numStepsInWindow, const scalar_t ts , aux_states_container_t & auxStatesContainer) const{
  int n = 0;
  scalar_t t = ts + n*dt;

  hess.setZero(); // may be possible to get rid of these
  gradient.setZero();

  const auto wlsCurrentState = ::pressio::containers::span(wlsState,0,romSize);
  fomStateReconstrObj_(wlsCurrentState,yFOM_current);
  updateStatesFirstStep(wlsStateIC,fomStateReconstrObj_,auxStatesContainer,ode);

  residualPolicy(ode,appObj,yFOM_current,residual,auxStatesContainer,ts,dt);
  for (int i = 0; i < time_stencil_size; i++){
    jacobianPolicy(ode,appObj,yFOM_current,wlsJacs[time_stencil_size - i - 1],phi_,auxStatesContainer,n*dt,dt,i);
  }
  hess_type C(romSize,romSize); // this is a temporary work around for the += issue.
  auto hess_block = ::pressio::containers::subspan( hess,std::make_pair( n*romSize,(n+1)*romSize ) , std::make_pair( n*romSize,(n+1)*romSize ) );
  ::pressio::containers::ops::dot(wlsJacs[time_stencil_size-1],wlsJacs[time_stencil_size-1], hess_block);
  local_vec_update(wlsJacs[time_stencil_size-1],residual,gradient,n*romSize,romSize);
  for (int n = 1; n < numStepsInWindow; n++){
    // === reconstruct FOM states ========
    updateStatesNStep(yFOM_current,auxStatesContainer,ode); 

    const auto wlsCurrentState = ::pressio::containers::span(wlsState,n*romSize,romSize);
    fomStateReconstrObj_(wlsCurrentState,yFOM_current);
    // == Evaluate residual ============
    t = ts + n*dt;
    residualPolicy(ode,appObj,yFOM_current,residual,auxStatesContainer,t,dt);
    // == Evaluate Jacobian (Currently have Jacobian w.r.p to previous state hard coded outside of loop)
    jacobianPolicy(ode,appObj,yFOM_current,wlsJacs[time_stencil_size-1],phi_,auxStatesContainer,t,dt,0);
    // == Update everything
    int sbar = std::min(n,time_stencil_size);
    for (int i=0; i < sbar; i++){
      local_vec_update(wlsJacs[time_stencil_size-i-1],residual,gradient,(n-i)*romSize,romSize);}
    // == Assemble local component of global Hessian //
    for (int i=0; i < sbar; i++){
      for (int j=0; j <= i; j++){
         auto hess_block = ::pressio::containers::subspan(hess, std::make_pair( (n-i)*romSize,(n-i+1)*romSize ) , std::make_pair( (n-j)*romSize,(n-j+1)*romSize ) );
         ::pressio::containers::ops::dot(wlsJacs[time_stencil_size-i-1],wlsJacs[time_stencil_size-j-1], C);
         for (int k = 0; k<romSize; k++){
           for (int l = 0; l<romSize; l++){
             hess_block(k,l) = hess_block(k,l) + C(k,l); 
            }
         }
         auto hess_block2 = ::pressio::containers::subspan(hess, std::make_pair( (n-j)*romSize,(n-j+1)*romSize ) , std::make_pair( (n-i)*romSize,(n-i+1)*romSize ) );
         for (int k = 0; k<romSize; k++){
           for (int l = 0; l<romSize; l++){
             hess_block2(l,k) = hess_block(k,l); 
            }
         }
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

