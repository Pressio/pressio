namespace pressio{ namespace rom{ namespace wls{ namespace impl{
/* 
This header file contains policies for computing the time local Jacobians; i.e., 
the the Jacobian at the current time step with respect to 
the differet states in the (time) stencil
*/
struct local_jacobian_policy_velocityAPI{
  // Policies for local residuals and jacobians
  //BDF2
  template <typename fom_type,
            typename fom_state_type, 
            typename jac_type, 
            typename basis_type, 
            typename aux_states_container_type, 
            typename scalar_type>
  void operator()(const ::pressio::ode::implicitmethods::BDF2 & tag, 
                  const fom_type & appObj,
                  const fom_state_type & yFOM, 
                  jac_type & Jphi,
                  const basis_type & phi, 
                  const aux_states_container_type & auxStatesContainer, 
                  const scalar_type & t, 
                  const scalar_type & dt, 
                  const int & step, 
                  int arg=0 ) const{
    // u^n - u^{n-1} - dt*f ;
    if (step == 0){ 
      if (arg == 0){
        appObj.applyJacobian(*yFOM.data(),*phi.data(),t,*(Jphi).data());
        constexpr auto cn   = ::pressio::ode::constants::bdf2<scalar_type>::c_n_; //      1
        ::pressio::containers::ops::do_update(Jphi,-dt,phi,cn);}  
      }
    if (step > 0){ 
      if (arg == 0){
        appObj.applyJacobian(*yFOM.data(),*phi.data(),t,*(Jphi).data());
        constexpr auto cn   = ::pressio::ode::constants::bdf2<scalar_type>::c_n_; //      1
        const auto cfdt   = ::pressio::ode::constants::bdf2<scalar_type>::c_f_*dt; //     2/3
        ::pressio::containers::ops::do_update(Jphi,cfdt,phi,cn);}
      if (arg == 1 && step == 0){//only perform computation once since this never changes
        constexpr auto cnm1   = ::pressio::ode::constants::bdf2<scalar_type>::c_nm1_; // -4/3
        ::pressio::containers::ops::do_update(Jphi,phi,cnm1);}
      if (arg == 2 && step == 0){//only perform computation once since this never changes
        constexpr auto cnm2   = ::pressio::ode::constants::bdf2<scalar_type>::c_nm2_; //  2/3
        ::pressio::containers::ops::do_update(Jphi,phi,cnm2);}
      }
  }
  //Implicit Euler
  template <typename fom_type, 
            typename fom_state_type, 
            typename jac_type, 
            typename basis_type, 
            typename aux_states_container_type, 
            typename scalar_type>
  void operator()(const ::pressio::ode::implicitmethods::Euler & tag,
                  const fom_type & appObj,
                  const fom_state_type & yFOM, 
                  jac_type & Jphi, 
                  const basis_type & phi, 
                  const aux_states_container_type & auxStatesContainer,
                  const scalar_type & t, 
                  const scalar_type & dt, 
                  const int & step, 
                  int arg=0 ) const{
    // u^n - u^{n-1} - f ;
    if (arg == 0){
      appObj.applyJacobian(*yFOM.data(),*phi.data(),t,*(Jphi).data());
      constexpr auto cn   = ::pressio::ode::constants::bdf1<scalar_type>::c_n_; //      1
      const auto cfdt   = ::pressio::ode::constants::bdf1<scalar_type>::c_f_*dt; //  -1*dt
      ::pressio::containers::ops::do_update(Jphi,cfdt,phi,cn);}
    if (arg == 1 && step == 0){//only perform computation once since this never changes
      constexpr auto cnm1   = ::pressio::ode::constants::bdf1<scalar_type>::c_nm1_; // -1.
      ::pressio::containers::ops::do_update(Jphi,phi,cnm1);}
  }

  //Explicit Euler
  template <typename fom_type, 
            typename fom_state_type, 
            typename jac_type, 
            typename basis_type, 
            typename aux_states_container_type, 
            typename scalar_type>
  void operator()(const ::pressio::ode::explicitmethods::Euler & tag,
                  const fom_type & appObj,
                  const fom_state_type & yFOM, 
                  jac_type & Jphi, 
                  const basis_type & phi, 
                  const aux_states_container_type & auxStatesContainer,
                  const scalar_type & t, 
                  const scalar_type & dt, 
                  const int & step, 
                  int arg=0 ) const{
    // u^n - u^{n-1} - f ;
    if (arg == 0 && step == 0){//only perform computation once since this never changes
      const auto One  = ::pressio::utils::constants::one<scalar_type>(); //  1*dt
      ::pressio::containers::ops::do_update(Jphi,phi,One);}
    if (arg == 1){
      using nm1 = ode::nMinusOne;
      const auto dtnegOne  = ::pressio::utils::constants::negOne<scalar_type>()*dt; //  -1*dt
      const auto One  = ::pressio::utils::constants::one<scalar_type>(); //  1*dt

      auto & odeState_nm1 = auxStatesContainer.get(nm1());
      appObj.applyJacobian(*odeState_nm1.data(),*phi.data(),t,*(Jphi).data());
      ::pressio::containers::ops::do_update(Jphi,dtnegOne,phi,One);}
  }


};
//============================
}}}}
