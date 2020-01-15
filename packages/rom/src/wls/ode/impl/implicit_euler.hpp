
namespace pressio{ namespace rom{ namespace wls{ namespace ode{ namespace impl{

template<typename fom_state_t>
class ImplicitEuler{
private:
  int & stateSize_;
  using nm1 = ::pressio::ode::nMinusOne;
public:

  ImplicitEuler(int & stateSize) : stateSize_(stateSize){}

  // Residual Policy
  template <typename fom_type,
          typename fom_state_type,
          typename residual_type,
          typename aux_states_container_type,
          typename scalar_type>
  void time_discrete_residual(
                  const fom_type & appObj,
                  const fom_state_type & yFOM,
                  residual_type & residual,
                  const aux_states_container_type & auxStatesContainer,
                  const scalar_type & t,
                  const scalar_type & dt,
                  const int & step) const{
    appObj.velocity(*yFOM.data(),t,*residual.data());
    ::pressio::ode::impl::time_discrete_residual<::pressio::ode::implicitmethods::Euler>(yFOM,residual,auxStatesContainer,dt);
  }

  // Jacobian Policy
  template <typename fom_type,
            typename fom_state_type,
            typename jac_type,
            typename basis_type,
            typename aux_states_container_type,
            typename scalar_type>
  void time_discrete_jacobian(
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

  // Policies to update containers

  // Policy for first step, here we move data from the IC container to the aux container
  template <typename wls_state_type,
            typename fom_state_reconstr_t,
            typename aux_states_container_t>
  void updateStatesFirstStep(const wls_state_type & wlsStateIC,
                             const fom_state_reconstr_t & fomStateReconstr,
                             aux_states_container_t & auxStatesContainer) const{
      const auto wlsInitialStateNm1 = ::pressio::containers::span(wlsStateIC,0,this->stateSize_);
      auto & fomStateNm1 = auxStatesContainer.get(nm1());
      fomStateReconstr(wlsInitialStateNm1,fomStateNm1);
  }

  // Policies at an N step. Here we are just dealing with the aux container
  template <typename aux_states_container_t>
  void updateStatesNStep(const fom_state_t & yFOM_current_,
                         aux_states_container_t & auxStatesContainer) const{
      auto & odeState_nm1 = auxStatesContainer.get(nm1());
      ::pressio::containers::ops::deep_copy(yFOM_current_, odeState_nm1);
  }


};

}}}}}
