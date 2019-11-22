namespace pressio{ namespace rom{ namespace wls{ namespace impl{

//struct JTJ_JTR_policy_standard{
//public:
//  void operator()(wls_state_type wls_state, fom_t appObj) const { 
//    std::cout << "Computed JTJ and JTR through default policy" << std::endl;
//}
//};

template< typename wls_state_type, typename fom_type, typename jtj_type, typename jtr_type, typename decoder_t>
struct JTJ_JTR_policy_smart{
public:
  using jtj_t = jtj_type;
  using jtr_t = jtr_type;
  using fom_native_state_t      = typename fom_type::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using fom_state_reconstr_t    = pressio::rom::FomStateReconstructor<fom_state_t, decoder_t>;
  //using fom_state_container_t = typename pressio::ode::StatesContainer<fom_state_t, 2>;
  using scalar_t	= typename fom_type::scalar_type;

  template <typename fom_state_container_t>
  void operator()(fom_type & appObj, wls_state_type  & wlsState, jtj_type & jtj, jtr_type & jtr,fom_state_container_t &  fomStateContainerObj_, fom_state_reconstr_t & fomStateReconstrObj_,scalar_t dt, std::size_t numStepsInWindow) const { 
  for (int i = 1; i < numStepsInWindow; i++){
    fomStateReconstrObj_(wlsState[i],fomStateContainerObj_[0]);
    fomStateReconstrObj_(wlsState[i-1],fomStateContainerObj_[1]);
    std::cout << fomStateContainerObj_[0][i]<< std::endl;
    appObj.timeDiscreteResidual(i,dt*i,dt,*fomStateContainerObj_[0].data(),*fomStateContainerObj_[0].data(),*fomStateContainerObj_[1].data());
    //appObj.applyTimeDiscreteJacobian(i,dt*i,dt,*fomStateContainerObj_[0].data(),*fomStateContainerObj_[0].data(),*fomStateContainerObj_[1].data());
    std::cout << "Computed JTJ and JTR through the smart policy" << std::endl;
  } 
  }
};



template< typename wls_state_type, typename fom_type, typename residual_type, typename decoder_t>
struct residual_policy_naive{
public:
  using residual_t = residual_type;
  using fom_native_state_t      = typename fom_type::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using fom_state_reconstr_t    = pressio::rom::FomStateReconstructor<fom_state_t, decoder_t>;
  using fom_state_container_t = typename pressio::ode::StatesContainer<fom_state_t, 2>;
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




}}}}

