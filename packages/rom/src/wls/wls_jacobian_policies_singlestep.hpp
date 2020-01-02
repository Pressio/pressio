namespace pressio{ namespace rom{ namespace wls{ namespace impl{
/* 
This header file contains policies for computing the time local Jacobians; i.e., 
the the Jacobian at the current time step with respect to 
the differet states in the (time) stencil
*/
struct local_jacobian_policy_velocityAPI{
public:
  // Policies for local residuals and jacobians
  //BDF2
  template <typename fom_type, 
						typename fom_state_t, 
						typename jac_t, 
						typename basis_t, 
						typename aux_states_container_t, 
						typename scalar_t>
  void operator()(const ::pressio::ode::implicitmethods::BDF2 & tag, 
									const fom_type & appObj,
									const fom_state_t & yFOM, 
									jac_t & Jphi,
									const basis_t & phi, 
									const aux_states_container_t & auxStatesContainer, 
									const scalar_t & t, 
									const scalar_t & dt, 
									const int & step, 
									int arg=0 ) const{
    // u^n - u^{n-1} - dt*f ;
    if (step == 0){ 
      if (arg == 0){
        appObj.applyJacobian(*yFOM.data(),*phi.data(),t,*(Jphi).data());
        ::pressio::containers::ops::do_update(Jphi,-dt,phi,1.);}  //::pressio::ode::impl::time_discrete_jacobian
        //::pressio::ode::impl::time_discrete_jacobian<::pressio::ode::implicitmethods::Euler>(Jphi,dt);}
      }
    if (step > 0){ 
      if (arg == 0){
        appObj.applyJacobian(*yFOM.data(),*phi.data(),t,*(Jphi).data());
        ::pressio::containers::ops::do_update(Jphi,-2./3.*dt,phi,1.);}
      if (arg == 1){
        ::pressio::containers::ops::do_update(Jphi,phi,-4./3.);}
      if (arg == 2){
        ::pressio::containers::ops::do_update(Jphi,phi,1./3.);}
      }
  }
  //Implicit Euler
  template <typename fom_type, 
						typename fom_state_t, 
						typename jac_t, 
						typename basis_t, 
						typename aux_states_container_t, 
						typename scalar_t>
  void operator()(const ::pressio::ode::implicitmethods::Euler & tag,
									const fom_type & appObj,
									const fom_state_t & yFOM, 
									jac_t & Jphi, 
									const basis_t & phi, 
									const aux_states_container_t & auxStatesContainer,
									const scalar_t & t, 
									const scalar_t & dt, 
									const int & step, 
									int arg=0 ) const{
    // u^n - u^{n-1} - f ;
    if (arg == 0){
      appObj.applyJacobian(*yFOM.data(),*phi.data(),t,*(Jphi).data());
      ::pressio::containers::ops::do_update(Jphi,-dt,phi,1.);}
    if (arg == 1){
      ::pressio::containers::ops::do_update(Jphi,phi,-1.);}
  }
};
//============================
}}}}
