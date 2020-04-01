/*
//@HEADER
// ************************************************************************
//
//  rom_wls_bdf2.hpprom_wls_bdf2.hpp
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

#ifndef ROM_WLS_BDF2_HPP_
#define ROM_WLS_BDF2_HPP_

namespace pressio{ namespace rom{ namespace wls{ namespace timeschemes{ namespace impl{

template<typename fom_state_t, typename wls_state_t>
class BDF2{

public:
  static constexpr int state_stencil_size_ = 2;
  static constexpr bool is_explicit	   = false;

private:
  using aux_states_container_t = ::pressio::ode::AuxStatesContainer<is_explicit, fom_state_t, state_stencil_size_>;

  ::pressio::rom::wls::rom_size_t romStateSize_;
  mutable aux_states_container_t auxStatesContainer_;
  mutable bool jacobianZeroNeedsRecomputing_ = true;
  mutable bool jacobianOneNeedsRecomputing_  = true;

public:
  BDF2() = delete;
  BDF2(const BDF2 &) = delete;
  BDF2(BDF2 &&) = delete;
  BDF2 & operator=(const BDF2 &) = delete;
  BDF2 & operator=(BDF2 &&) = delete;

  BDF2(int & romStateSize,const fom_state_t & fomState)
    :  romStateSize_(romStateSize),
       auxStatesContainer_(fomState)
  {}

public:
  template <
    typename fom_type,
    typename fom_state_type,
    typename residual_type,
    typename scalar_type>
  void time_discrete_residual(const fom_type & appObj,
			      const fom_state_type & fomState,
			      residual_type & residual,
			      const scalar_type & t,
			      const  scalar_type & dt,
			      const int & step) const
  {
    if (step > 0){
      // u^n - 4./3.*u^{n-1} + 1./3.u^{n-2} - 2./3.*dt*f
      appObj.velocity(*fomState.data(), t, *residual.data());
      ::pressio::ode::impl::time_discrete_residual<
	::pressio::ode::implicitmethods::BDF2>(fomState,residual,auxStatesContainer_,dt);
    }

    if (step == 0){
      // u^n - u^{n-1}  - dt*f
      appObj.velocity(*fomState.data(), t, *residual.data());
      ::pressio::ode::impl::time_discrete_residual<
	::pressio::ode::implicitmethods::Euler>(fomState, residual, auxStatesContainer_, dt);
    }
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
			      const int & step,
			      int arg=0 ) const
  {

    // u^n - u^{n-1} - dt*f ;
    if (step == 0){
      if (arg == 0){
        appObj.applyJacobian(*fomState.data(), *phi.data(), t, *(Jphi).data());
        constexpr auto cn   = ::pressio::ode::constants::bdf2<scalar_type>::c_n_; //1
        ::pressio::ops::do_update(Jphi, -dt, phi, cn);
      }
    }

    if (step > 0){
      if (arg == 0){
        appObj.applyJacobian(*fomState.data(),*phi.data(),t,*(Jphi).data());
        constexpr auto cn   = ::pressio::ode::constants::bdf2<scalar_type>::c_n_; // 1
        const auto cfdt   = ::pressio::ode::constants::bdf2<scalar_type>::c_f_*dt; //2/3
        ::pressio::ops::do_update(Jphi, cfdt, phi, cn);
      }

      if (arg == 1 && jacobianOneNeedsRecomputing_){//only perform computation once since this never changes
        constexpr auto cnm1   = ::pressio::ode::constants::bdf2<scalar_type>::c_nm1_; // -4/3
        ::pressio::ops::do_update(Jphi, phi, cnm1);
        jacobianOneNeedsRecomputing_ = false;
        std::cout << "here" << std::endl;
      }

      if (arg == 2 && jacobianZeroNeedsRecomputing_){//only perform computation once since this never changes
        constexpr auto cnm2   = ::pressio::ode::constants::bdf2<scalar_type>::c_nm2_; //  2/3
        ::pressio::ops::do_update(Jphi, phi, cnm2);
        jacobianZeroNeedsRecomputing_ = false;
        std::cout << "here" << std::endl;
      }
    }
  }

  bool jacobianNeedsRecomputing(std::size_t i) const{
    return (i==0) ? jacobianZeroNeedsRecomputing_ : jacobianOneNeedsRecomputing_;
  }


  // for first step, here we move data from the IC container to the aux container
  template <typename wls_state_type,typename fom_state_reconstr_t>
  void updateStatesFirstStep(const wls_state_type & wlsStateIC,
                             const fom_state_reconstr_t & fomStateReconstr) const
  {
    using nm1 = ::pressio::ode::nMinusOne;
    using nm2 = ::pressio::ode::nMinusTwo;

    const auto wlsInitialStateNm2 = ::pressio::containers::span(wlsStateIC, 0, romStateSize_);
    const auto wlsInitialStateNm1 = ::pressio::containers::span(wlsStateIC, romStateSize_, romStateSize_);
    auto & fomStateNm2 = auxStatesContainer_.get(nm2());
    auto & fomStateNm1 = auxStatesContainer_.get(nm1());
    fomStateReconstr(wlsInitialStateNm2,fomStateNm2);
    fomStateReconstr(wlsInitialStateNm1,fomStateNm1);
  }

  // at an N step. Here we are just dealing with the aux container
  void updateStatesNStep(const fom_state_t & fomStateCurrent) const
  {
    using nm1 = ::pressio::ode::nMinusOne;
    using nm2 = ::pressio::ode::nMinusTwo;

    auto & odeState_nm1 = auxStatesContainer_.get(nm1());
    auto & odeState_nm2 = auxStatesContainer_.get(nm2());
    ::pressio::ops::deep_copy(odeState_nm2, odeState_nm1);
    ::pressio::ops::deep_copy(odeState_nm1, fomStateCurrent);
  }

};

}}}}} // end namespace pressio::rom::wls::ode::impl
#endif
