namespace pressio{ namespace rom{ namespace wls{ namespace impl{

//struct JTJ_JTR_policy_standard{
//public:
//  void operator()(wls_state_type wls_state, fom_t appObj) const { 
//    std::cout << "Computed JTJ and JTR through default policy" << std::endl;
//}
//};

template< typename wls_state_type, typename fom_type, typename hess_type, typename gradient_type, typename decoder_t>
struct JTJ_JTR_policy_smart{
public:
  using hess_t = hess_type;
  using gradient_t = gradient_type;
  using fom_native_state_t      = typename fom_type::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using fom_state_reconstr_t    = pressio::rom::FomStateReconstructor<fom_state_t, decoder_t>;
  //using fom_state_container_t = typename pressio::ode::StatesContainer<fom_state_t, 2>;
  using residual_t = fom_state_t;
  using scalar_t	= typename fom_type::scalar_type;

  using jac_t  = typename decoder_t::jacobian_t;
  using wls_jacs_t     =  std::vector<jac_t>;

  // Fixed Types
  using eig_dyn_mat = Eigen::Matrix<scalar_t, -1, -1>;
  using jtj_t = pressio::containers::MultiVector<eig_dyn_mat>;
  
  jac_t & phi_;
  mutable wls_jacs_t wlsJacs;
  hess_t JTJ;
  hess_t JTJ1;
  mutable residual_t residual;
  int romSize; 
  int time_stencil_size;
  int fomSize;
  // Policy constructor
  JTJ_JTR_policy_smart(int time_stencil_size,jac_t & phi, int romSize, int numStepsInWindow, int fomSize) : wlsJacs(time_stencil_size,phi), phi_(phi), JTJ(romSize*numStepsInWindow,romSize*numStepsInWindow),JTJ1(romSize,romSize), residual(fomSize)  {
    this->romSize = romSize;    
    this->time_stencil_size = time_stencil_size;
    this->fomSize = fomSize;
  }

  template <typename fom_state_container_t>
  void operator()(fom_type & appObj, wls_state_type  & wlsState, wls_state_type & wlsStateIC, hess_type & jtj, gradient_type & gradient,fom_state_container_t &  fomStateContainerObj_, fom_state_reconstr_t & fomStateReconstrObj_,scalar_t dt, std::size_t numStepsInWindow) const { 
  int n = 0;


//  const auto wlsCurrentState = wlsState.viewColumnVector(n)
  fomStateReconstrObj_(wlsState[n],fomStateContainerObj_[1]);





  fomStateReconstrObj_(wlsStateIC[0],fomStateContainerObj_[0]);
  (*wlsJacs[0].data()).block(0,0,fomSize,romSize) = -(*phi_.data());
  appObj.timeDiscreteResidual(n,dt*n,dt,*residual.data(),*fomStateContainerObj_[1].data(),*fomStateContainerObj_[0].data());
  appObj.applyTimeDiscreteJacobian(n,dt*n,dt,*phi_.data(),0,*(wlsJacs[1]).data(),*fomStateContainerObj_[1].data(), *fomStateContainerObj_[0].data());
  (*jtj.data()).block(n*romSize,n*romSize,romSize,romSize) = (*wlsJacs[1].data()).transpose() * (*wlsJacs[1].data());


  (*gradient.data()).block(n*romSize,0,romSize,1) = -(*wlsJacs[1].data()).transpose()*(*residual.data());


  for (int n = 1; n < numStepsInWindow; n++){
    // === reconstruct FOM states ========
    *fomStateContainerObj_[0].data() = *fomStateContainerObj_[1].data();
    fomStateReconstrObj_(wlsState[n],fomStateContainerObj_[1]);
    // == Evaluate residual ============
    appObj.timeDiscreteResidual(n,dt*n,dt,*residual.data(),*fomStateContainerObj_[1].data(),*fomStateContainerObj_[0].data());
    // == Evaluate Jacobian (Currently have jacobian w.r.p to previous state hard coded outside of loop)
    appObj.applyTimeDiscreteJacobian(n,dt*n,dt,*phi_.data(),0,*(wlsJacs[1]).data(),*fomStateContainerObj_[1].data(), *fomStateContainerObj_[0].data());
    (*gradient.data()).block(n*romSize,0,romSize,1) = -(*wlsJacs[1].data()).transpose()*(*residual.data());
    for (int i=1; i < time_stencil_size; i++){
      (*gradient.data()).block((n-i)*romSize,0,romSize,1) -= (*wlsJacs[time_stencil_size-i-1].data()).transpose()*(*residual.data());}
    // == Assembel local component of global Hessian //
    for (int i=0; i < time_stencil_size; i++){
      for (int j=0; j <= i; j++){
         (*jtj.data()).block((n-i)*romSize,(n-j)*romSize,romSize,romSize) += (*wlsJacs[time_stencil_size-i-1].data()).transpose() * (*wlsJacs[time_stencil_size-j-1].data());
         if (i != j){
          (*jtj.data()).block((n-j)*romSize,(n-i)*romSize,romSize,romSize) =  (*jtj.data()).block((n-i)*romSize,(n-j)*romSize,romSize,romSize).transpose();}
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

