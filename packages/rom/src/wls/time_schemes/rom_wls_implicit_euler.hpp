/*
//@HEADER
// ************************************************************************
//
// rom_wls_implicit_euler.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef ROM_WLS_TIME_SCHEMES_ROM_WLS_IMPLICIT_EULER_HPP_
#define ROM_WLS_TIME_SCHEMES_ROM_WLS_IMPLICIT_EULER_HPP_

namespace pressio{ namespace rom{ namespace wls{ namespace timeschemes{ namespace impl{

template<typename fom_state_t>
class ImplicitEuler
{

public:
  static constexpr ::pressio::rom::wls::window_size_t state_stencil_size_ = 1;

private:
  ::pressio::rom::wls::rom_size_t romStateSize_;
  using aux_states_container_t = ::pressio::ode::implicitmethods::StencilStatesManager
    <fom_state_t, state_stencil_size_>;
  mutable aux_states_container_t stencilStates_;
  mutable bool jacobianNeedsRecomputing_ = true;

public:
  ImplicitEuler() = delete;
  ImplicitEuler(const ImplicitEuler &) = delete;
  ImplicitEuler(ImplicitEuler &&) = delete;
  ImplicitEuler & operator=(const ImplicitEuler &) = delete;
  ImplicitEuler & operator=(ImplicitEuler &&) = delete;

  ImplicitEuler(::pressio::rom::wls::rom_size_t & romStateSize,
		const fom_state_t & fomState)
    : romStateSize_(romStateSize),
      stencilStates_(fomState)
  {}

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  ImplicitEuler(::pressio::rom::wls::rom_size_t & romStateSize,
		const typename fom_state_t::traits::wrapped_t & fomState)
    : romStateSize_(romStateSize),
      stencilStates_(fom_state_t(fomState))
  {}
#endif

public:
  bool jacobianNeedsRecomputing(std::size_t i) const{
    return jacobianNeedsRecomputing_;
  }

  // for first step: move data from the IC container to the aux container
  template <typename wls_state_type, typename fom_state_reconstr_t>
  void updateStatesFirstStep(const wls_state_type & wlsStateIC,
                             const fom_state_reconstr_t & fomStateReconstr) const
  {
    const auto wlsInitialState = ::pressio::containers::span(wlsStateIC, 0, this->romStateSize_);
    auto & fomState = stencilStates_(::pressio::ode::n());
    fomStateReconstr(wlsInitialState,fomState);
  }

  // at an N step we are just dealing with the aux container
  void updateStatesNStep(const fom_state_t & fomStateCurrent) const
  {
    auto & odeState = stencilStates_(::pressio::ode::n());
    ::pressio::ops::deep_copy(odeState, fomStateCurrent);
  }

  /* ==================================================================
			time-continuous API
     ================================================================== */

#ifdef PRESSIO_ENABLE_TPL_EIGEN
  /** time_discrete_residual for eigen */
  template <
    typename fom_type,
    typename fom_state_type,
    typename residual_type,
    typename scalar_type
    >
  ::pressio::mpl::enable_if_t<
    ::pressio::rom::constraints::most_likely_continuous_time_system<fom_type>::value and
    ::pressio::containers::predicates::is_vector_wrapper_eigen<fom_state_type>::value == true
    >
  time_discrete_residual(const fom_type & fomSystemObj,
			 const fom_state_type & fomState,
			 residual_type & residual,
			 const scalar_type & t,
			 const scalar_type & dt,
			 const window_size_t & step) const
  {
    fomSystemObj.velocity(*fomState.data(),t,*residual.data());
    ::pressio::ode::impl::discrete_time_residual(fomState, residual,
						 stencilStates_, dt,
						 ::pressio::ode::implicitmethods::Euler());
  }

  /** time_discrete_jacobian for eigen */
  template <
    typename fom_type,
    typename fom_state_type,
    typename jac_type,
    typename basis_type,
    typename scalar_type
    >
  ::pressio::mpl::enable_if_t<
    ::pressio::rom::constraints::most_likely_continuous_time_system<fom_type>::value and
    ::pressio::containers::predicates::is_vector_wrapper_eigen<fom_state_type>::value == true
    >
  time_discrete_jacobian(const fom_type & fomSystemObj,
			 const fom_state_type & fomState,
			 jac_type & Jphi,
			 const basis_type & phi,
			 const scalar_type & t,
			 const scalar_type & dt,
			 const ::pressio::rom::wls::window_size_t & step,
			 int arg=0 ) const
  {
    // u^n - u^{n-1} - f ;
    if (arg == 0){
      fomSystemObj.applyJacobian(*fomState.data(), *phi.data(), t, *(Jphi).data());
      constexpr auto cnp1   = ::pressio::ode::constants::bdf1<scalar_type>::c_np1_;
      const auto cfdt     = ::pressio::ode::constants::bdf1<scalar_type>::c_f_*dt; //  -1*dt
      ::pressio::ops::update(Jphi, cfdt, phi, cnp1);
    }

    //only perform computation once since this never changes
    if (arg == 1 && jacobianNeedsRecomputing_){
      constexpr auto cn   = ::pressio::ode::constants::bdf1<scalar_type>::c_n_; // -1.
      ::pressio::ops::update(Jphi, phi, cn);
      jacobianNeedsRecomputing_ = true;
    }
  }
#endif

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  template <
    typename fom_type,
    typename fom_state_type,
    typename residual_type,
    typename scalar_type
    >
  ::pressio::mpl::enable_if_t<
    ::pressio::rom::constraints::most_likely_continuous_time_system<fom_type>::value and
    ::pressio::containers::predicates::is_rank1_tensor_wrapper_pybind<fom_state_type>::value == true
    >
  time_discrete_residual(const fom_type & fomSystemObj,
			 const fom_state_type & fomState,
			 residual_type & residual,
			 const scalar_type & t,
			 const scalar_type & dt,
			 const window_size_t & step) const
  {
    fomSystemObj.velocity(*fomState.data(),t,*residual.data());
    ::pressio::ode::impl::discrete_time_residual(fomState, residual,
						 stencilStates_, dt,
						 ::pressio::ode::implicitmethods::Euler());
  }

  template <
    typename fom_type,
    typename fom_state_type,
    typename jac_type,
    typename basis_type,
    typename scalar_type
    >
  ::pressio::mpl::enable_if_t<
    ::pressio::rom::constraints::most_likely_continuous_time_system<fom_type>::value and
    ::pressio::containers::predicates::is_rank1_tensor_wrapper_pybind<fom_state_type>::value == true
    >
  time_discrete_jacobian(const fom_type & fomSystemObj,
			 const fom_state_type & fomState,
			 jac_type & Jphi,
			 const basis_type & phi,
			 const scalar_type & t,
			 const scalar_type & dt,
			 const ::pressio::rom::wls::window_size_t & step,
			 int arg=0 ) const
  {
    // u^n - u^{n-1} - f ;
    if (arg == 0){
      fomSystemObj.applyJacobian(*fomState.data(), *phi.data(), t, *(Jphi).data());
      constexpr auto cnp1   = ::pressio::ode::constants::bdf1<scalar_type>::c_np1_;
      const auto cfdt     = ::pressio::ode::constants::bdf1<scalar_type>::c_f_*dt; //  -1*dt
      ::pressio::ops::update(Jphi, cfdt, phi, cnp1);
    }

    //only perform computation once since this never changes
    if (arg == 1 && jacobianNeedsRecomputing_){
      constexpr auto cn   = ::pressio::ode::constants::bdf1<scalar_type>::c_n_; // -1.
      ::pressio::ops::update(Jphi, phi, cn);
      jacobianNeedsRecomputing_ = true;
    }
  }
#endif


#ifdef PRESSIO_ENABLE_TPL_KOKKOS
  /** time_discrete_residual for kokkos */
  template <
    typename fom_type,
    typename fom_state_type,
    typename residual_type,
    typename scalar_type
    >
  ::pressio::mpl::enable_if_t<
    ::pressio::rom::constraints::most_likely_continuous_time_system<fom_type>::value and
    ::pressio::containers::predicates::is_vector_wrapper_kokkos<fom_state_type>::value == true
    >
  time_discrete_residual(const fom_type & fomSystemObj,
			 const fom_state_type & fomState,
			 residual_type & residual,
			 const scalar_type & t,
			 const scalar_type & dt,
			 const window_size_t & step) const
  {
    fomSystemObj.velocity(*fomState.data(),t,*residual.data());
    ::pressio::ode::impl::discrete_time_residual(fomState, residual,
						 stencilStates_, dt,
						 ::pressio::ode::implicitmethods::Euler());
  }

  /** time_discrete_jacobian for kokkos */
  template <
    typename fom_type,
    typename fom_state_type,
    typename jac_type,
    typename basis_type,
    typename scalar_type
    >
  ::pressio::mpl::enable_if_t<
    ::pressio::rom::constraints::most_likely_continuous_time_system<fom_type>::value and
    ::pressio::containers::predicates::is_vector_wrapper_kokkos<fom_state_type>::value == true
    >
  time_discrete_jacobian(const fom_type & fomSystemObj,
			 const fom_state_type & fomState,
			 jac_type & Jphi,
			 const basis_type & phi,
			 const scalar_type & t,
			 const scalar_type & dt,
			 const ::pressio::rom::wls::window_size_t & step,
			 int arg=0 ) const
  {
    // u^n - u^{n-1} - f ;
    if (arg == 0){
      fomSystemObj.applyJacobian(*fomState.data(), *phi.data(), t, *(Jphi).data());
      constexpr auto cnp1   = ::pressio::ode::constants::bdf1<scalar_type>::c_np1_; //      1
      const auto cfdt     = ::pressio::ode::constants::bdf1<scalar_type>::c_f_*dt; //  -1*dt
      ::pressio::ops::update(Jphi, cfdt, phi, cnp1);
    }

    //only perform computation once since this never changes
    if (arg == 1 && jacobianNeedsRecomputing_){
      constexpr auto cn   = ::pressio::ode::constants::bdf1<scalar_type>::c_n_; // -1.
      ::pressio::ops::update(Jphi, phi, cn);
      jacobianNeedsRecomputing_ = true;
    }
  }
#endif

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
  /* time-discrete residual w/ tpetra */
  template <
    typename fom_type,
    typename fom_state_type,
    typename residual_type,
    typename scalar_type
    >
  ::pressio::mpl::enable_if_t<
    ::pressio::rom::constraints::most_likely_continuous_time_system<fom_type>::value and
    ::pressio::containers::predicates::is_vector_wrapper_tpetra<fom_state_type>::value==true
    >
  time_discrete_residual(const fom_type & fomSystemObj,
			 const fom_state_type & fomState,
			 residual_type & residual,
			 const scalar_type & t,
			 const scalar_type & dt,
			 const window_size_t & step) const
  {
    fomSystemObj.velocity(*fomState.data(),t,*residual.data());
    auto residualView = *residual.data();
    auto fomStateView = *fomState.data();
    auto Un = stencilStates_(::pressio::ode::n());
    auto UnView = *Un.data();
    _time_discrete_residual_from_views(fomStateView,UnView,residualView,step,dt);
  }

  /* time-discrete residual tpetra block */
  template <
    typename fom_type,
    typename fom_state_type,
    typename residual_type,
    typename scalar_type
    >
  ::pressio::mpl::enable_if_t<
    ::pressio::rom::constraints::most_likely_continuous_time_system<fom_type>::value and
    ::pressio::containers::predicates::is_vector_wrapper_tpetra_block<fom_state_type>::value == true
    >
  time_discrete_residual(const fom_type & fomSystemObj,
			 const fom_state_type & fomState,
			 residual_type & residual,
			 const scalar_type & t,
			 const scalar_type & dt,
			 const window_size_t & step) const
  {
    fomSystemObj.velocity(*fomState.data(),t,*residual.data());
    auto residualNative = *residual.data();
    auto residualView = residualNative.getVectorView();
    auto fomStateNative = *fomState.data();
    auto fomStateView = fomStateNative.getVectorView();
    auto Un = stencilStates_(::pressio::ode::n());
    auto UnView = (*Un.data()).getVectorView();
    _time_discrete_residual_from_views(fomStateView,UnView,residualView,step,dt);
  }

  /* time_discrete_jacobian w/ tpetra block */
  template <
    typename fom_type,
    typename fom_state_type,
    typename jac_type,
    typename basis_type,
    typename scalar_type
    >
  ::pressio::mpl::enable_if_t<
    ::pressio::rom::constraints::most_likely_continuous_time_system<fom_type>::value and
    ::pressio::containers::predicates::is_vector_wrapper_tpetra_block<fom_state_type>::value == true
    >
  time_discrete_jacobian(const fom_type & fomSystemObj,
			 const fom_state_type & fomState,
			 jac_type & Jphi,
			 const basis_type & phi,
			 const scalar_type & t,
			 const scalar_type & dt,
			 const ::pressio::rom::wls::window_size_t & step,
			 int arg=0 ) const
  {

    auto JphiNative = *Jphi.data();
    auto JphiView = JphiNative.getMultiVectorView();
    auto phiNative = *phi.data();
    auto phiView = phiNative.getMultiVectorView();
    auto fomStateNative = *fomState.data();
    auto fomStateView = fomStateNative.getVectorView();

    _time_discrete_jacobian_from_views(fomSystemObj,fomState,fomStateView,Jphi,
				       JphiView,phi,phiView,t,dt,step,arg);
  }

  /* time_discrete_jacobian w/ tpetra */
  template <
    typename fom_type,
    typename fom_state_type,
    typename jac_type,
    typename basis_type,
    typename scalar_type
    >
  ::pressio::mpl::enable_if_t<
    ::pressio::rom::constraints::most_likely_continuous_time_system<fom_type>::value and
    ::pressio::containers::predicates::is_vector_wrapper_tpetra<fom_state_type>::value == true
    >
  time_discrete_jacobian(const fom_type & fomSystemObj,
			 const fom_state_type & fomState,
			 jac_type & Jphi,
			 const basis_type & phi,
			 const scalar_type & t,
			 const scalar_type & dt,
			 const ::pressio::rom::wls::window_size_t & step,
			 int arg=0 ) const
  {
    auto JphiView = *Jphi.data();
    auto phiView = *phi.data();
    auto fomStateView = *fomState.data();
    _time_discrete_jacobian_from_views(fomSystemObj,fomState,fomStateView,Jphi,
				       JphiView,phi,phiView,t,dt,step,arg);
  }
#endif


  /* ==================================================================
			time-discrete API
     ================================================================== */

  /* time_discrete_residual */
  template <
    typename fom_type,
    typename fom_state_type,
    typename residual_type,
    typename scalar_type
    >
  ::pressio::mpl::enable_if_t<
    ::pressio::rom::constraints::most_likely_discrete_time_system<fom_type>::value
    >
  time_discrete_residual(const fom_type & fomSystemObj,
			 const fom_state_type & fomState,
			 residual_type & residual,
			 const scalar_type & t,
			 const scalar_type & dt,
			 const window_size_t & step) const
  {
    auto & fomStateN = stencilStates_(::pressio::ode::n());
    fomSystemObj.discreteTimeResidual(step, t, dt,
				      *residual.data(),
				      *fomState.data(), *fomStateN.data());
  }

  /* time_discrete_jacobian */
  template <
    typename fom_type,
    typename fom_state_type,
    typename jac_type,
    typename basis_type,
    typename scalar_type
    >
  ::pressio::mpl::enable_if_t<
    ::pressio::rom::constraints::most_likely_discrete_time_system<fom_type>::value
    >
  time_discrete_jacobian(const fom_type & fomSystemObj,
			 const fom_state_type & fomState,
			 jac_type & Jphi,
			 const basis_type & phi,
			 const scalar_type & t,
			 const scalar_type & dt,
			 const ::pressio::rom::wls::window_size_t & step,
			 int arg=0 ) const
  {
    // u^n+1 - u^{n} - f ;
    if (arg == 0){
      auto & fomStateN = stencilStates_(::pressio::ode::n());
      fomSystemObj.applyDiscreteTimeJacobian(step, t, dt, *phi.data(),
					     *Jphi.data(), *fomState.data(),
					     *fomStateN.data());
    }

    //only perform computation once since this never changes
    if (arg == 1 && jacobianNeedsRecomputing_){
      constexpr auto cn   = ::pressio::ode::constants::bdf1<scalar_type>::c_n_; // -1.
      ::pressio::ops::update(Jphi,phi,cn);
      jacobianNeedsRecomputing_ = true;
    }
  }

private:
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
  // time scheem loop for tpetra(block) data structures
  template<
  typename scalar_type,
  typename fom_state_view_t,
  typename residualView_t,
  typename step_t
  >
  void _time_discrete_residual_from_views(const fom_state_view_t & fomStateView,
					 const fom_state_view_t & Unm1View,
					 residualView_t & residualView,
					 const step_t & step,
					 const scalar_type & dt) const
  {
    const auto hyperMap = residualView.getMap();
    const auto gIDr = hyperMap->getMyGlobalIndices();
    const auto fomStateMap = fomStateView.getMap();
    const auto gIDy = fomStateMap->getMyGlobalIndices();
    auto Rdv = residualView.getDataNonConst();
    auto Udv = fomStateView.getData();
    auto Unm1dv = Unm1View.getData();
    const auto cfdt     = ::pressio::ode::constants::bdf1<scalar_type>::c_f_*dt; //  -1*dt
    // get my global elements
    for (size_t i=0; i<residualView.getLocalLength(); i++){
      const auto lid = fomStateMap->getLocalElement(gIDr[i]);
      auto Rd = Rdv[i];
      auto Ud = Udv[lid];
      auto Unm1d = Unm1dv[lid];
      Rd = Ud - Unm1d + cfdt*Rd;
      residualView.replaceLocalValue(i,Rd);
    }
  }


  /**
     time_discrete_jacobian function overload for continuous time API w/ tpetra block
  */
  template <
    typename fom_type,
    typename fom_state_type,
    typename jac_type,
    typename basis_type,
    typename fom_state_view_type,
    typename jac_view_type,
    typename basis_view_type,
    typename scalar_type
    >
  ::pressio::mpl::enable_if_t<
    ::pressio::rom::constraints::most_likely_continuous_time_system<fom_type>::value
    >
  _time_discrete_jacobian_from_views(const fom_type & fomSystemObj,
				     const fom_state_type & fomState,
				     const fom_state_view_type & fomStateView,
				     jac_type & Jphi,
				     jac_view_type & JphiView,
				     const basis_type & phi,
				     const basis_view_type & phiView,
				     const scalar_type & t,
				     const scalar_type & dt,
				     const ::pressio::rom::wls::window_size_t & step,
				     int arg=0 ) const
  {
    const auto hyperMap = JphiView.getMap();
    const auto gIDJphi = hyperMap->getMyGlobalIndices();
    const auto fomStateMap = fomStateView.getMap();
    const auto gIDy = fomStateMap->getMyGlobalIndices();

    if (arg == 0){
      fomSystemObj.applyJacobian(*fomState.data(), *phi.data(), t, *(Jphi).data());
      const auto cfdt     = ::pressio::ode::constants::bdf1<scalar_type>::c_f_*dt; //  -1*dt

      // get my global elements
      for (size_t i=0; i<JphiView.getLocalLength(); i++){
        const auto lid = fomStateMap->getLocalElement(gIDJphi[i]);
        for (size_t k=0 ; k < (size_t)Jphi.extent(1); k++){
          auto Jphid = JphiView.getDataNonConst(k)[i];
          auto phid = phiView.getData(k)[lid];
          Jphid = cfdt*Jphid + phid;
          JphiView.replaceLocalValue(i,k,Jphid);
        }
      }
    }

    //only perform computation once since this never changes
    if (arg == 1 && jacobianNeedsRecomputing_){
      constexpr auto cn   = ::pressio::ode::constants::bdf1<scalar_type>::c_n_; // -1.
      // get my global elements
      for (size_t i=0; i<JphiView.getLocalLength(); i++){
        for (size_t k=0 ; k < (size_t)Jphi.extent(1); k++){
          const auto lid = fomStateMap->getLocalElement(gIDJphi[i]);
          auto phid = phiView.getData(k)[lid];
          JphiView.replaceLocalValue(i,k,cn*phid);
        }
      }
      jacobianNeedsRecomputing_ = true;
    }
  }
#endif
};

}}}}} // end namespace pressio::rom::wls::ode::impl
#endif  // ROM_WLS_TIME_SCHEMES_ROM_WLS_IMPLICIT_EULER_HPP_
