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

#ifndef ROM_WLS_IMPLICIT_EULER_HPP_
#define ROM_WLS_IMPLICIT_EULER_HPP_

namespace pressio{ namespace rom{ namespace wls{ namespace timeschemes{ namespace impl{

template<typename fom_state_t, typename wls_state_t>
class ImplicitEuler{

public:
  static constexpr ::pressio::rom::wls::window_size_t state_stencil_size_ = 1;
  static constexpr bool is_explicit	   = false;

private:
  ::pressio::rom::wls::rom_size_t romStateSize_;
  using aux_states_container_t = ::pressio::ode::AuxStatesContainer<is_explicit, fom_state_t, state_stencil_size_>;
  mutable aux_states_container_t auxStatesContainer_;
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
      auxStatesContainer_(fomState)
  {}

public:
  template <
    typename fom_type,
    typename fom_state_type,
    typename residual_type,
    typename scalar_type
    >
  void time_discrete_residual(const fom_type & appObj,
			      const fom_state_type & fomState,
			      residual_type & residual,
			      const scalar_type & t,
			      const scalar_type & dt,
			      const window_size_t & step) const
  {
    appObj.velocity(*fomState.data(),t,*residual.data());
    ::pressio::ode::impl::time_discrete_residual<
      ::pressio::ode::implicitmethods::Euler
      >(fomState, residual, auxStatesContainer_, dt);
  }

  template <
    typename fom_type,
    typename fom_state_type,
    typename jac_type,
    typename basis_type,
    typename scalar_type>
  void time_discrete_jacobian(const fom_type & appObj,
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
      appObj.applyJacobian(*fomState.data(),*phi.data(),t,*(Jphi).data());
      constexpr auto cn   = ::pressio::ode::constants::bdf1<scalar_type>::c_n_; //      1
      const auto cfdt     = ::pressio::ode::constants::bdf1<scalar_type>::c_f_*dt; //  -1*dt
      ::pressio::ops::do_update(Jphi,cfdt,phi,cn);
    }

    //only perform computation once since this never changes
    if (arg == 1 && jacobianNeedsRecomputing_){
      constexpr auto cnm1   = ::pressio::ode::constants::bdf1<scalar_type>::c_nm1_; // -1.
      ::pressio::ops::do_update(Jphi,phi,cnm1);
      jacobianNeedsRecomputing_ = false;
    }
  }

  bool jacobianNeedsRecomputing(std::size_t i) const{
    return jacobianNeedsRecomputing_;
  }

  // for first step: move data from the IC container to the aux container
  template <typename wls_state_type, typename fom_state_reconstr_t>
  void updateStatesFirstStep(const wls_state_type & wlsStateIC,
                             const fom_state_reconstr_t & fomStateReconstr) const
  {
    using nm1 = ::pressio::ode::nMinusOne;
    const auto wlsInitialStateNm1 = ::pressio::containers::span(wlsStateIC, 0, this->romStateSize_);
    auto & fomStateNm1 = auxStatesContainer_.get(nm1());
    fomStateReconstr(wlsInitialStateNm1,fomStateNm1);
  }

  // at an N step we are just dealing with the aux container
  void updateStatesNStep(const fom_state_t & fomStateCurrent) const
  {
    using nm1 = ::pressio::ode::nMinusOne;
    auto & odeState_nm1 = auxStatesContainer_.get(nm1());
    ::pressio::ops::deep_copy(odeState_nm1, fomStateCurrent);
  }

};

}}}}} // end namespace pressio::rom::wls::ode::impl
#endif
