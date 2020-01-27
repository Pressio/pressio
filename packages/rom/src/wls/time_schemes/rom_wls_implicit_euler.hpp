
#ifndef ROM_WLS_IMPLICIT_EULER_HPP_
#define ROM_WLS_IMPLICIT_EULER_HPP_

namespace pressio{ namespace rom{ namespace wls{ namespace timeschemes{ namespace impl{

template<typename fom_state_t, typename wls_state_t>
class ImplicitEuler{

public:
  static constexpr int state_stencil_size_ = 2;
  static constexpr bool is_explicit	   = false;

private:
  int stateSize_;
  using aux_states_container_t = ::pressio::ode::AuxStatesContainer<is_explicit, fom_state_t, state_stencil_size_>;
  mutable aux_states_container_t auxStatesContainer_;

public:
  ImplicitEuler() = delete;
  ImplicitEuler(const ImplicitEuler &) = delete;
  ImplicitEuler(ImplicitEuler &&) = delete;
  ImplicitEuler & operator=(const ImplicitEuler &) = delete;
  ImplicitEuler & operator=(ImplicitEuler &&) = delete;

  ImplicitEuler(int & stateSize, const fom_state_t & yFOM)
    : stateSize_(stateSize),
      auxStatesContainer_(yFOM)
  {}

  template <
    typename fom_type,
    typename fom_state_type,
    typename residual_type,
    typename scalar_type
    >
  void time_discrete_residual(const fom_type & appObj,
			      const fom_state_type & yFOM,
			      residual_type & residual,
			      const scalar_type & t,
			      const scalar_type & dt,
			      const int & step) const
  {
    appObj.velocity(*yFOM.data(),t,*residual.data());
    ::pressio::ode::impl::time_discrete_residual<
      ::pressio::ode::implicitmethods::Euler>(yFOM, residual, auxStatesContainer_, dt);
  }

  template <
    typename fom_type,
    typename fom_state_type,
    typename jac_type,
    typename basis_type,
    typename scalar_type>
  void time_discrete_jacobian(const fom_type & appObj,
			      const fom_state_type & yFOM,
			      jac_type & Jphi,
			      const basis_type & phi,
			      const scalar_type & t,
			      const scalar_type & dt,
			      const int & step,
			      int arg=0 ) const
  {
    // u^n - u^{n-1} - f ;
    if (arg == 0){
      appObj.applyJacobian(*yFOM.data(),*phi.data(),t,*(Jphi).data());
      constexpr auto cn   = ::pressio::ode::constants::bdf1<scalar_type>::c_n_; //      1
      const auto cfdt     = ::pressio::ode::constants::bdf1<scalar_type>::c_f_*dt; //  -1*dt
      ::pressio::containers::ops::do_update(Jphi,cfdt,phi,cn);
    }
    if (arg == 1 && step == 0){//only perform computation once since this never changes
      constexpr auto cnm1   = ::pressio::ode::constants::bdf1<scalar_type>::c_nm1_; // -1.
      ::pressio::containers::ops::do_update(Jphi,phi,cnm1);
    }
  }

  // for first step: move data from the IC container to the aux container
  template <typename wls_state_type, typename fom_state_reconstr_t>
  void updateStatesFirstStep(const wls_state_type & wlsStateIC,
                             const fom_state_reconstr_t & fomStateReconstr) const
  {
    using nm1 = ::pressio::ode::nMinusOne;
    const auto wlsInitialStateNm1 = ::pressio::containers::span(wlsStateIC,0,this->stateSize_);
    auto & fomStateNm1 = auxStatesContainer_.get(nm1());
    fomStateReconstr(wlsInitialStateNm1,fomStateNm1);
  }

  // at an N step we are just dealing with the aux container
  void updateStatesNStep(const fom_state_t & yFOM_current_) const
  {
    using nm1 = ::pressio::ode::nMinusOne;
    auto & odeState_nm1 = auxStatesContainer_.get(nm1());
    ::pressio::containers::ops::deep_copy(yFOM_current_, odeState_nm1);
  }

};

}}}}} // end namespace pressio::rom::wls::ode::impl
#endif
