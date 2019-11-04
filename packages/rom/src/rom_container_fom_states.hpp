/*
//@HEADER
// ************************************************************************
//
// rom_container_fom_states.hpp
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

#ifndef ROM_CONTAINER_FOM_STATES_HPP_
#define ROM_CONTAINER_FOM_STATES_HPP_

#include "rom_ConfigDefs.hpp"
#include "rom_fwd.hpp"
#include "../../ode/src/ode_states_container.hpp"

namespace pressio{ namespace rom{

template <
  typename fom_state_type,
  int N,
  typename reconstuctor_type
  >
class FomStatesContainer<fom_state_type, N, reconstuctor_type>
{
public:
  static constexpr int N_ = N;

  // fom_state_type & operator[](std::size_t i){
  //   assert( i<N );
  //   fomOldStates_[i];
  // }

  // fom_state_type const & operator[](std::size_t i) const{
  //   assert( i<N );
  //   fomOldStates_[i];
  // }

  FomStatesContainer() = delete;
  // FomStatesContainer(const FomStatesContainer &) = delete;
  // FomStatesContainer & operator=(const FomStatesContainer &) = delete;
  // FomStatesContainer(FomStatesContainer &&) = delete;
  // FomStatesContainer & operator=(FomStatesContainer &&) = delete;
  ~FomStatesContainer() = default;

  /* ----------------
   * N = 0
   * ---------------*/

  /* cnstr for N = 0, and fom_state_type is a pressio vector wrapper */
  template <
    typename _fom_state_type = fom_state_type,
    int _N = N,
    ::pressio::mpl::enable_if_t<
      _N==0 and
      ::pressio::containers::meta::is_vector_wrapper<_fom_state_type>::value
      > * = nullptr
    >
  FomStatesContainer(const _fom_state_type & fomStateIn,
		const reconstuctor_type & fomStateReconstr)
    : fomState_(fomStateIn),
      fomStateReconstrObj_(fomStateReconstr)
  {
    this->resetContainersToZero();
  }

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  /* cnstr for N = 0, and fom_state_type is a pybind11 array */
  template <
    typename _fom_state_type = fom_state_type,
    int _N = N,
    ::pressio::mpl::enable_if_t<
      _N==0 and
      ::pressio::containers::meta::is_array_pybind11<_fom_state_type>::value
      > * = nullptr
    >
  FomStatesContainer(const _fom_state_type & fomStateIn,
		const reconstuctor_type & fomStateReconstr)
    : fomState_{ {_fom_state_type(const_cast<_fom_state_type &>(fomStateIn).request())} },
      fomStateReconstrObj_(fomStateReconstr)
  {
    this->resetContainersToZero();
  }
#endif


  /* ----------------
   * N = 1
   * ---------------*/
  /* cnstr for N = 1, and fom_state_type is a pressio vector wrapper */
  template <
    typename _fom_state_type = fom_state_type,
    int _N = N,
    ::pressio::mpl::enable_if_t<
      _N==1 and
      ::pressio::containers::meta::is_vector_wrapper<_fom_state_type>::value
      > * = nullptr
    >
  FomStatesContainer(const _fom_state_type & fomStateIn,
		const reconstuctor_type & fomStateReconstr)
    : fomState_(fomStateIn),
      fomOldStates_{{fomStateIn}},
      fomStateReconstrObj_(fomStateReconstr)
  {
    this->resetContainersToZero();
  }

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  /* cnstr for N = 1, and fom_state_type is a pybind11 array */
  template <
    typename _fom_state_type = fom_state_type,
    int _N = N,
    ::pressio::mpl::enable_if_t<
      _N==1 and
      ::pressio::containers::meta::is_array_pybind11<_fom_state_type>::value
      > * = nullptr
    >
  FomStatesContainer(const _fom_state_type & fomStateIn,
		const reconstuctor_type & fomStateReconstr)
    : fomState_{ {_fom_state_type(const_cast<_fom_state_type &>(fomStateIn).request())} },
      fomOldStates_{ {_fom_state_type(const_cast<_fom_state_type &>(fomStateIn).request())} },
      fomStateReconstrObj_(fomStateReconstr)
  {
    this->resetContainersToZero();
  }
#endif


  /* ----------------
   * N = 2
   * ---------------*/
  /* cnstr for N = 2, and fom_state_type is a pressio vector wrapper */
  template <
    typename _fom_state_type = fom_state_type,
    int _N = N,
    ::pressio::mpl::enable_if_t<
      _N==2 and
      ::pressio::containers::meta::is_vector_wrapper<_fom_state_type>::value
      > * = nullptr
    >
  FomStatesContainer(const _fom_state_type & fomStateIn,
		const reconstuctor_type & fomStateReconstr)
    : fomState_(fomStateIn),
      fomOldStates_{{fomStateIn, fomStateIn}},
      fomStateReconstrObj_(fomStateReconstr)
  {
    this->resetContainersToZero();
  }

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  /* cnstr for N = 2, and fom_state_type is a pybind11 array */
  template <
    typename _fom_state_type = fom_state_type,
    int _N = N,
    ::pressio::mpl::enable_if_t<
      _N==2 and
      ::pressio::containers::meta::is_array_pybind11<_fom_state_type>::value
      > * = nullptr
    >
  FomStatesContainer(const _fom_state_type & fomStateIn,
		const reconstuctor_type & fomStateReconstr)
    : fomState_{ {_fom_state_type(const_cast<_fom_state_type &>(fomStateIn).request())} },
      fomOldStates_{ {_fom_state_type(const_cast<_fom_state_type &>(fomStateIn).request())},
		{_fom_state_type(const_cast<_fom_state_type &>(fomStateIn).request())}},
      fomStateReconstrObj_(fomStateReconstr)
  {
    this->resetContainersToZero();
  }
#endif


public:

  const fom_state_type & getCRefToCurrentFomState() const{
    return fomState_;
  }

  const std::array<fom_state_type, N> & getCRefToFomOldStates() const{
    return fomOldStates_;
  }

  template <typename rom_state_t>
  void reconstructCurrentFomState(const rom_state_t & romY) const
  {
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("reconstruct fom state");
#endif

    fomStateReconstrObj_(romY, fomState_);

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("reconstruct fom state");
#endif
  }

  template <int n, typename rom_state_t>
  void reconstructFomOldStates(const ::pressio::ode::StatesContainer<rom_state_t, n> & romYprev) const
  {
#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("reconstruct fom old state");
#endif

    for (auto i=0; i<n; i++){
      fomStateReconstrObj_(romYprev[i], fomOldStates_[i]);
    }

#ifdef PRESSIO_ENABLE_TEUCHOS_TIMERS
    timer->stop("reconstruct fom old state");
#endif
  }

private:
  /* set all entries to zero for all members */
  void resetContainersToZero(){
    ::pressio::containers::ops::set_zero(fomState_);
    for (auto i=0; i<N; i++)
      ::pressio::containers::ops::set_zero(fomOldStates_[i]);
  }

private:
  mutable fom_state_type fomState_                = {};
  mutable std::array<fom_state_type, N> fomOldStates_  = {};
  const reconstuctor_type & fomStateReconstrObj_  = {};

};//end class

}}//end namespace pressio::rom
#endif
