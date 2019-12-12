namespace pressio{ namespace rom{ namespace wls{ namespace impl{
/* 
This header file contains policies for computing the time local Jacobians and residuals; i.e., 
the residual at a current time step, and the Jacobian at the current time step with respect to 
the differet states in the (time) stencil
*/

//============ For the velocity APIs ===========================
template< typename fom_type>
struct local_residual_policy_velocityAPI{
public:
  using fom_native_state_t      = typename fom_type::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using fom_state_container_t = std::vector<fom_state_t>;
  using scalar_t	= typename fom_type::scalar_type;
  using residual_t = fom_state_t; //maybe we need to change this for hyper reduction?
  void operator()(fom_type & appObj, residual_t & residual,fom_state_container_t & fomStateContainerObj_,scalar_t t, scalar_t dt) const{
    // u^n - u^{n-1} - f
    appObj.velocity(*fomStateContainerObj_[1].data(),t,*residual.data());
    (*residual.data()) = *(fomStateContainerObj_[1]).data() - *(fomStateContainerObj_[0]).data() - dt*(*residual.data());
    }
};

template< typename fom_type>
struct local_residual_policy_BDF2_velocityAPI{
public:
  using fom_native_state_t      = typename fom_type::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using fom_state_container_t = std::vector<fom_state_t>;
  using scalar_t	= typename fom_type::scalar_type;
  using residual_t = fom_state_t; //maybe we need to change this for hyper reduction?
  void operator()(fom_type & appObj, residual_t & residual,fom_state_container_t & fomStateContainerObj_,scalar_t t, scalar_t dt) const{
    if (t > 1e-10){
    // u^n - 4./3.*u^{n-1} + 1./3.u^{n-1} - 2./3.*dt*f
    appObj.velocity(*fomStateContainerObj_[2].data(),t,*residual.data());
    (*residual.data()) = *(fomStateContainerObj_[2]).data() - 4./3.* *(fomStateContainerObj_[1]).data() + 1./3.* *(fomStateContainerObj_[0]).data()  - 2./3.*dt*(*residual.data());} 
    if (t < 1e-10){
    // u^n - 4./3.*u^{n-1} + 1./3.u^{n-1} - 2./3.*dt*f
    std::cout << "Here" << std::endl;
    appObj.velocity(*fomStateContainerObj_[2].data(),t,*residual.data());
    (*residual.data()) = *(fomStateContainerObj_[2]).data() - *(fomStateContainerObj_[1]).data()  - dt*(*residual.data());} 
    }
};


template< typename fom_type, typename decoder_t>
struct local_jacobian_policy_velocityAPI{
public:
  using fom_native_state_t      = typename fom_type::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using fom_state_container_t = std::vector<fom_state_t>;
  using scalar_t	= typename fom_type::scalar_type;
  using basis_t  = typename decoder_t::jacobian_t;
  using jac_t = typename pressio::containers::Matrix<Eigen::Matrix<scalar_t, -1, -1>>;
  // Policies for local residuals and jacobians
  void operator()(fom_type & appObj, jac_t & Jphi, basis_t & phi, fom_state_container_t & fomStateContainerObj_,scalar_t t, scalar_t dt, int arg=0) const{
    // u^n - u^{n-1} - f ; 
    if (arg == 0){
      appObj.applyJacobian(*fomStateContainerObj_[1].data(),*phi.data(),t,*(Jphi).data());
      (*Jphi.data()) = *phi.data() - dt*(*Jphi.data());}
    if (arg == 1){
      (*Jphi.data()) = -*phi.data() ;}
    }
};


template< typename fom_type, typename decoder_t>
struct local_jacobian_policy_BDF2_velocityAPI{
public:
  using fom_native_state_t      = typename fom_type::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using fom_state_container_t = std::vector<fom_state_t>;
  using scalar_t	= typename fom_type::scalar_type;
  using basis_t  = typename decoder_t::jacobian_t;
  using jac_t = typename pressio::containers::Matrix<Eigen::Matrix<scalar_t, -1, -1>>;
  // Policies for local residuals and jacobians
  void operator()(fom_type & appObj, jac_t & Jphi, basis_t & phi, fom_state_container_t & fomStateContainerObj_,scalar_t t, scalar_t dt, int arg=0) const{
    // u^n - u^{n-1} - f ;
    if (t < 1e-10){ 
      if (arg == 0){
        appObj.applyJacobian(*fomStateContainerObj_[1].data(),*phi.data(),t,*(Jphi).data());
        (*Jphi.data()) = *phi.data() - dt*(*Jphi.data());}
      if (arg == 1){
        (*Jphi.data()) = -4./3.*(*phi.data()) ;}
      if (arg == 2){
        (*Jphi.data()) =  1./3.*(*phi.data()) ;}
      }
    if (t > 1e-10){ 
      if (arg == 0){
        appObj.applyJacobian(*fomStateContainerObj_[1].data(),*phi.data(),t,*(Jphi).data());
        (*Jphi.data()) = *phi.data() - 2./3.*dt*(*Jphi.data());}
      if (arg == 1){
        (*Jphi.data()) = -4./3.*(*phi.data()) ;}
      if (arg == 2){
        (*Jphi.data()) =  1./3.*(*phi.data()) ;}
      }

  }
};


// ===================================================================

//============ For the residual APIs ====================================
template< typename fom_type>
struct local_residual_policy_residualAPI{
public:
  using fom_native_state_t      = typename fom_type::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using fom_state_container_t = std::vector<fom_state_t>;
  using scalar_t	= typename fom_type::scalar_type;
  using residual_t = fom_state_t; //maybe we need to change this for hyper reduction?
  void operator()(fom_type & appObj, residual_t & residual,fom_state_container_t & fomStateContainerObj_,scalar_t t, scalar_t dt) const{
    // u^n - u^{n-1} - f
    int n = (int) t/dt;
    appObj.timeDiscreteResidual(n,t,dt,*residual.data(),*fomStateContainerObj_[1].data(),*fomStateContainerObj_[0].data());
    }
};

template< typename fom_type, typename decoder_t>
struct local_jacobian_policy_residualAPI{
public:
  using fom_native_state_t      = typename fom_type::state_type;
  using fom_state_t             = ::pressio::containers::Vector<fom_native_state_t>;
  using fom_state_container_t = std::vector<fom_state_t>;
  using scalar_t	= typename fom_type::scalar_type;
  using basis_t  = typename decoder_t::jacobian_t;
  using jac_t = typename pressio::containers::Matrix<Eigen::Matrix<scalar_t, -1, -1>>;
  // Policies for local residuals and jacobians
  void operator()(fom_type & appObj, jac_t & Jphi, basis_t & phi, fom_state_container_t & fomStateContainerObj_,scalar_t t, scalar_t dt, int arg=0) const{
    // u^n - u^{n-1} - f ; 
    if (arg == 0){
      int n = (int) t/dt;
      appObj.applyTimeDiscreteJacobian(n,dt*n,dt,*phi.data(),0,*Jphi.data(),*fomStateContainerObj_[1].data(), *fomStateContainerObj_[0].data());
    }
    if (arg == 1){
      (*Jphi.data()) = -*phi.data() ;}
    }
};
//============================
}}}}
