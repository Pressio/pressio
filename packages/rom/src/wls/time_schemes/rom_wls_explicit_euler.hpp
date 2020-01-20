
#ifndef ROM_WLS_EXPLICIT_EULER_HPP_
#define ROM_WLS_EXPLICIT_EULER_HPP_

namespace pressio{ namespace rom{ namespace wls{ namespace timeschemes{ namespace impl{

template<typename fom_state_t, typename wls_state_t>
class ExplicitEuler{

public:
  static constexpr int state_stencil_size_ = 2;
  static constexpr bool is_explicit	   = true;

private:
  using aux_states_container_t = ::pressio::ode::AuxStatesContainer<is_explicit, fom_state_t, state_stencil_size_>;
  int stateSize_;
  mutable aux_states_container_t auxStatesContainer_;

public:
  ExplicitEuler() = delete;
  ExplicitEuler(const ExplicitEuler &) = delete;
  ExplicitEuler(ExplicitEuler &&) = delete;
  ExplicitEuler & operator=(const ExplicitEuler &) = delete;
  ExplicitEuler & operator=(ExplicitEuler &&) = delete;

  ExplicitEuler(int & stateSize, const fom_state_t & yFOM)
    : stateSize_(stateSize),
      auxStatesContainer_(yFOM)
  {}

  template <
    typename fom_type,
    typename fom_state_type,
    typename residual_type,
    typename scalar_type>
  void time_discrete_residual(const fom_type & appObj,
			      const fom_state_type & yFOM,
			      residual_type & residual,
			      const scalar_type & t,
			      const scalar_type & dt,
			      const int & step) const
  {
    // auto & odeState_nm1 = auxStatesContainer.get(nm1());
    // const auto dtnegOne  = ::pressio::utils::constants::negOne<scalar_type>()*dt; //  -1*dt
    // constexpr auto negOne  = ::pressio::utils::constants::negOne<scalar_type>(); //  -1*dt
    // constexpr auto  One  = ::pressio::utils::constants::one<scalar_type>(); //  -1*dt
    // appObj.velocity(*odeState_nm1.data(),t,*residual.data()); // gives f
    // ::pressio::containers::ops::do_update(residual,dtnegOne,yFOM,One,odeState_nm1,negOne);
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
    // // u^n - u^{n-1} - f ;
    // if (arg == 0 && step == 0){//only perform computation once since this never changes
    //   const auto One  = ::pressio::utils::constants::one<scalar_type>(); //  1*dt
    //   ::pressio::containers::ops::do_update(Jphi,phi,One);
    // }

    // if (arg == 1){
    //   const auto dtnegOne  = ::pressio::utils::constants::negOne<scalar_type>()*dt; //  -1*dt
    //   const auto One  = ::pressio::utils::constants::one<scalar_type>(); //  1*dt
    //   auto & odeState_nm1 = auxStatesContainer.get(nm1());
    //   appObj.applyJacobian(*odeState_nm1.data(),*phi.data(),t,*(Jphi).data());
    //   ::pressio::containers::ops::do_update(Jphi,dtnegOne,phi,One);
    // }
  }

  // first step, here we move data from the IC container to the aux container
  template <typename wls_state_type, typename fom_state_reconstr_t>
  void updateStatesFirstStep(const wls_state_type & wlsStateIC,
                             const fom_state_reconstr_t & fomStateReconstr) const
  {
    // const auto wlsInitialStateNm1 = ::pressio::containers::span(wlsStateIC,0,this->stateSize_);
    // auto & fomStateNm1 = auxStatesContainer_.get(nm1());
    // fomStateReconstr(wlsInitialStateNm1,fomStateNm1);
  }
  // at an N step. Here we are just dealing with the aux container
  void updateStatesNStep(const fom_state_t & yFOM_current_) const
  {
    // auto & odeState_nm1 = auxStatesContainer_.get(nm1());
    // ::pressio::containers::ops::deep_copy(yFOM_current_, odeState_nm1);
  }

};

}}}}} // end namespace pressio::rom::wls::ode::impl
#endif
